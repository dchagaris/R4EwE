
# Silence codetools NOTE about internal function fn.objfxn2
if (getRversion() >= "3.4.0") {
  utils::globalVariables(c("fn.objfxn2"))
}
