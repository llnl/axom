These classes instantiate ``axom::mir::ClipField`` for available devices and mesh
types. They were moved to the ``future`` directory for now since they can take
significant time to compile.

======================
ClipFieldFilterDevice
======================

This class instantiates ClipField for relevant mesh/coordset types on a specific device.

======================
ClipFieldFilter
======================

This class instantiates ClipFieldFilterDevice for all supported device types and
makes the policy selectable at runtime.

