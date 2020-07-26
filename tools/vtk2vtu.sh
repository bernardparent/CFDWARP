#!/bin/sh
#
# This is free and unencumbered software released into the public domain.
#
# Anyone is free to copy, modify, publish, use, compile, sell, or
# distribute this software, either in source code form or as a compiled
# binary, for any purpose, commercial or non-commercial, and by any
# means.
#
# In jurisdictions that recognize copyright laws, the author or authors
# of this software dedicate any and all copyright interest in the
# software to the public domain. We make this dedication for the benefit
# of the public at large and to the detriment of our heirs and
# successors. We intend this dedication to be an overt act of
# relinquishment in perpetuity of all present and future rights to this
# software under copyright law.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.
#
# For more information, please refer to <http://unlicense.org>


command -v pvpython >/dev/null 2>&1 || { echo "vtk2vtu.sh requires pvpython (part of the paraview distribution) but it is not installed.  Aborting." >&2; exit 1; }


if [ $# -gt 0 ]
then
{

echo "import sys
from paraview.simple import *

print 'Number of files to be converted:', len(sys.argv) - 1

for x in range(1, len(sys.argv)):        
    inputFile = str(sys.argv[x])
    outputFile = inputFile[:-1] + 'u'
    print x,': Converting ', inputFile, '  ->  ', outputFile
    reader = LegacyVTKReader( FileNames= inputFile )
    writer = XMLUnstructuredGridWriter()
    writer.FileName = outputFile
    writer.UpdatePipeline()
    Delete(reader)
    Delete(writer)" > .vtk2vtu.py


pvpython .vtk2vtu.py "$@"

rm -f .vtk2vtu.py

}
else
{
  echo ""
  echo "Converter of file format vtk to file format vtu using paraview python scripting pvpython."
  echo "Eg:"
  echo ""
  echo "$ vtk2vtu.sh post.vtk"
  echo "or"
  echo "$ vtk2vtu.sh post.vtk post2.vtk post3.vtk"
  echo ""
  echo "will create the post files post.vtu, post2.vtu, etc from the post files post.vtk, post2.vtk, etc."
  echo ""
} fi;

