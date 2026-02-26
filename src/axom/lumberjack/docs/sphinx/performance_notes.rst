Performance Notes: Tag-Selective Combining
=========================================

Context / Problem
-----------------

Some applications log very large numbers of messages with tags that are *not* intended
to be combined (e.g., millions of unique messages).  When Lumberjack is configured with
one or more combiners, :cpp:func:`axom::lumberjack::Lumberjack::combineMessages()` historically
performed an expensive nested scan across all messages to search for duplicates.

In the "many pass-through messages" case, this becomes the dominant cost even though
no combining should occur for those tags.

Goals and Constraints
---------------------

The desired behavior is:

* Allow different combining strategies for different tags.
* Preserve a single global sort-by-timestamp for messages written to a file/stream.
* Avoid paying quadratic duplicate-checking cost for tags that should never combine.
* Avoid multiplying MPI communication overhead by creating many Lumberjack instances
  (each Lumberjack implies its own reduction tree / communication pattern).

Design Thinking (Why This Approach)
-----------------------------------

The most direct way to avoid the quadratic scan is to ensure that messages that will never
be combined **do not participate** in the duplicate-check loops at all.

There are two broad ways to do that:

1. **Route messages into separate Lumberjack instances** based on tag ownership so that
   some Lumberjacks have no combiners and return early.  This can be effective but increases
   MPI cost because each Lumberjack performs its own tree reduction.
2. **Keep a single Lumberjack** and teach combiners to declare which messages they are
   willing to consider.  This keeps a single MPI tree while letting Lumberjack skip
   duplicate checking for pass-through messages.

The implemented change takes option (2) because it is the least invasive way to preserve
MPI behavior while addressing the observed bottleneck.

Implementation Summary
----------------------

Combiner hook
^^^^^^^^^^^^

`src/axom/lumberjack/Combiner.hpp` adds a new virtual hook:

* ``virtual bool isMessageCandidateForCombiner(const Message&)`` (default returns ``true``)

This allows a combiner to restrict itself to a subset of messages (commonly by tag).
If a message is not a candidate for **any** combiner, Lumberjack treats it as pass-through
and skips the expensive duplicate-checking loops for that message.

Updated combineMessages algorithm
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`src/axom/lumberjack/Lumberjack.cpp` updates ``combineMessages()`` to:

* Partition messages into:
  * **pass-through** messages (no combiner candidates) and
  * **candidate** messages (at least one combiner candidate)
* Run the nested duplicate-checking scan only for candidate messages.
* Append the candidate results back to the pass-through list.

Backwards compatibility:

* Existing combiners do not need to change because the default
  ``isMessageCandidateForCombiner()`` returns ``true``.
* With that default, behavior should match the prior algorithm (all messages are candidates),
  with only minor overhead from the extra candidate checks.

Guidance for tag-based combiners
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If a combiner uses tag ownership, implement:

* ``isMessageCandidateForCombiner(m)`` returning ``true`` only for owned tags.

Important requirement:

* If two messages could be combined by ``shouldMessagesBeCombined()``, then both must return
  ``true`` for ``isMessageCandidateForCombiner()``. Otherwise, a valid combine opportunity
  could be skipped.

How This Differs from ``feature/han12/lj_combiner_tags``
--------------------------------------------------------

This change was inspired by the same underlying idea (combiners should declare candidacy),
but differs in a few key ways:

* **No new Lumberjack tag state:** ``lj_combiner_tags`` introduced Lumberjack-managed tag lists
  (e.g., uncombinable tags) and explicit routing logic.  The implemented change does not add
  any Lumberjack-side configuration/state; candidacy is a combiner concern.
* **No per-combiner message bins:** ``lj_combiner_tags`` creates a separate ``finalMessages`` vector
  for each combiner and assigns each message to the first combiner that claims it.  The implemented
  change keeps one candidate list and one pass-through list, and allows a message to be a candidate
  for multiple combiners (it builds a candidate combiner list per message).
* **More direct "skip the scan" behavior:** The new algorithm explicitly short-circuits the nested
  duplicate-checking loops for messages that are not candidates for any combiner.

Targeted Perf Reproducer
------------------------

`src/axom/lumberjack/tests/lumberjack_speedTestTagCandidates.cpp` is a small standalone
test/benchmark that models the problematic pattern:

* A large number of pass-through messages with a tag that no combiner claims.
* A smaller number of messages with a tag that *is* combinable.

It is intended to make the performance issue and improvement obvious without needing to
set up multiple Lumberjack streams or MPI runs.

