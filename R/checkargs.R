# Internal TraMineR function for checking for deprecated arguments

checkargs <- function (arg.pairs) 
{
  new.names <- names(arg.pairs)
  old.names <- as.character(arg.pairs)
  pf <- parent.frame()
  fun.args <- ls(pf)
  for (i in 1:length(arg.pairs)) {
    new.name <- new.names[i]
    old.name <- old.names[i]
    if (new.name %in% fun.args && old.name %in% fun.args) {
      calling.fun <- paste0(as.character(sys.calls()[[1]])[1], 
                            "()")
      new.val <- eval(substitute(pf$nn, list(nn = new.name)))
      old.val <- eval(substitute(pf$on, list(on = old.name)))
      has.new <- !missing(new.val)
      has.old <- !missing(old.val)
      fun.args.default <- formals(sys.function(-1))
      new.default <- eval(substitute(fun.args.default$nn, 
                                     list(nn = new.name)))
      has.new.default <- !missing(new.default)
      if (has.new.default) {
        new.default <- tryCatch(eval(new.default), error = function(e) eval(substitute(pf$nd, 
                                                                                       list(nd = new.default))))
      }
      if (has.old) {
        if (has.new && (!has.new.default || (has.new.default && 
                                             !identical(new.val, new.default)))) {
          msg.stop("In", calling.fun, ":", new.name, 
                   "and", old.name, "cannot be specified together.", 
                   old.name, "is deprecated, use", new.name, 
                   "instead.")
        }
        else {
          msg.warn("In", calling.fun, ":", old.name, 
                   "is deprecated, use", new.name, "instead.")
          assign(new.name, old.val, envir = pf)
        }
      }
    }
    else {
      stop(new.name, " = ", old.name, " is incorrect: at least one argument doesn't exist")
    }
  }
}