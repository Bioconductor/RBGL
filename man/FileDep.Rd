\name{FileDep}
\alias{FileDep}
\non_function{}
\title{ FileDep: a graf object representing a file dependency dataset
example in boost graph library }
\usage{data(FileDep)}
\format{
  The format is:
 list()
 - attr(*, "vertices")= chr [1:15] "dax\_h" "yow\_h" "boz\_h" "zow\_h" ...
 - attr(*, "edges")=List of 19
  ..$ : chr [1:2] "dax\_h" "foo\_cpp"
  ..$ : chr [1:2] "dax\_h" "bar\_cpp"
  ..$ : chr [1:2] "dax\_h" "yow\_h"
  ..$ : chr [1:2] "yow\_h" "bar\_cpp"
  ..$ : chr [1:2] "yow\_h" "zag\_cpp"
  ..$ : chr [1:2] "boz\_h" "bar\_cpp"
  ..$ : chr [1:2] "boz\_h" "zig\_cpp"
  ..$ : chr [1:2] "boz\_h" "zag\_cpp"
  ..$ : chr [1:2] "zow\_h" "foo\_cpp"
  ..$ : chr [1:2] "foo\_cpp" "foo\_o"
  ..$ : chr [1:2] "foo\_o" "libfoobar\_a"
  ..$ : chr [1:2] "bar\_cpp" "bar\_o"
  ..$ : chr [1:2] "bar\_o" "libfoobar\_a"
  ..$ : chr [1:2] "libfoobar\_a" "libzigzag\_a"
  ..$ : chr [1:2] "zig\_cpp" "zig\_o"
  ..$ : chr [1:2] "zig\_o" "libzigzag\_a"
  ..$ : chr [1:2] "zag\_cpp" "zag\_o"
  ..$ : chr [1:2] "zag\_o" "libzigzag\_a"
  ..$ : chr [1:2] "libzigzag\_a" "killerapp"
  ..- attr(*, "costs")= symbol NULL
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
