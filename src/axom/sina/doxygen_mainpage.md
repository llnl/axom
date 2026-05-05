Sina {#sinatop}
=========

[Sina](@ref axom::sina) provides an easy way to collect data directly within codes and output them to a common file format. This is accomplished in an object oriented manner through the following classes:

- [Curve](@ref axom::sina::Curve): represents a 1D curve
- [CurveSet](@ref axom::sina::CurveSet): represents an entry in a record's "curve_set"
- [DataHolder](@ref axom::sina::DataHolder): a basic container for certain types of information
- [Datum](@ref axom::sina::Datum): tracks the value and (optionally) tags and/or units of a value associated with a Record
- [Document](@ref axom::sina::Document): represents the top-lvevl object of a JSON file conforming to the Sina schema
- [File](@ref axom::sina::File): tracks the location (URI) and mimetype of a file on the file system, plus any tags
- [Record](@ref axom::sina::Record): entry in a Document's Record list
- [Relationship](@ref axom::sina::Relationship): represents correlations between records; consists of three parts: a subject, an object, and a predicate
- [Run](@ref axom::sina::Run): a subtype of Record corresponding to a single run of an application, as specified in the Sina schema