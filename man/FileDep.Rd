\name{FileDep}
\alias{FileDep}
\non_function{}
\title{ FileDep: a graf object representing a file dependency dataset
example in boost graph library }
\usage{data(FileDep)}
\format{
  The format is:
 list()
 - attr(*, "vertices")= chr [1:15] "dax_h" "yow_h" "boz_h" "zow_h" ...
 - attr(*, "edges")=List of 19
  ..$ : chr [1:2] "dax_h" "foo_cpp"
  ..$ : chr [1:2] "dax_h" "bar_cpp"
  ..$ : chr [1:2] "dax_h" "yow_h"
  ..$ : chr [1:2] "yow_h" "bar_cpp"
  ..$ : chr [1:2] "yow_h" "zag_cpp"
  ..$ : chr [1:2] "boz_h" "bar_cpp"
  ..$ : chr [1:2] "boz_h" "zig_cpp"
  ..$ : chr [1:2] "boz_h" "zag_cpp"
  ..$ : chr [1:2] "zow_h" "foo_cpp"
  ..$ : chr [1:2] "foo_cpp" "foo_o"
  ..$ : chr [1:2] "foo_o" "libfoobar_a"
  ..$ : chr [1:2] "bar_cpp" "bar_o"
  ..$ : chr [1:2] "bar_o" "libfoobar_a"
  ..$ : chr [1:2] "libfoobar_a" "libzigzag_a"
  ..$ : chr [1:2] "zig_cpp" "zig_o"
  ..$ : chr [1:2] "zig_o" "libzigzag_a"
  ..$ : chr [1:2] "zag_cpp" "zag_o"
  ..$ : chr [1:2] "zag_o" "libzigzag_a"
  ..$ : chr [1:2] "libzigzag_a" "killerapp"
  ..- attr(*, "costs")= symbol NULL
  ..- attr(*, "class")= chr "edgeStruct"
 - attr(*, "accL")= list()
 - attr(*, "adjL")= list()
 - attr(*, "class")= chr "graf"
}
\source{
Boost Graph Library 
}
\references{
}
