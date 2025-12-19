Acknowledgement
================

Certain clipping tables have been borrowed/adapted from the VisIt software.

 - clipping/ClipCasesHex.cpp  
 - clipping/ClipCasesPyr.cpp
 - clipping/ClipCasesQua.cpp
 - clipping/ClipCasesTet.cpp
 - clipping/ClipCasesTri.cpp
 - clipping/ClipCasesWdg.cpp
 - split/ClipCasesPoly5.cpp
 - split/ClipCasesPoly6.cpp
 - split/ClipCasesPoly7.cpp
 - split/ClipCasesPoly8.cpp
 - split/ClipCasesQua.cpp

Regenerating Tables
====================

Polygonal clipping tables and cutting tables are derived from VisIt's
clipping files in the "split" directory. The name of the "split" directory
refers to the fact that the zone fragments are all split into triangles
and quads rather than appearing as polygons with more sides. The script
in split/convert_clip_cases.py can regenerate the polygonal clip and
cutting tables for polygon cases 4-8.

To regenerate the cutting tables:

	cd split
	python3 ./convert_clip_cases.py

To regenerate the clipping tables:

	cd split
	python3 ./convert_clip_cases.py -clip

Some light editing of the regenerated files may be required.


VisIt License
==============

The VisIt license appears below:

```
BSD 3-Clause License

Copyright (c) 2000 - 2025, Lawrence Livermore National Security, LLC
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```
