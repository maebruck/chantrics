# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

#' The \code{hello} function.
#'
#' Function printing \code{Hello, world!} in the console.
#'
#' Long long long long long long longlong long long long long long long long
#' long long long long long long long long long long
#'
#' @section Even more details:
#' Blah blah blah blah blah blah blah blah blah blah blah blah blah blah blah blah
#' @param times The number of times the phrase should be repeated.
#'
#' @return Returns nothing and prints \code{Hello, world!} to the console.
#' @seealso \code{\link[base]{rep}} combined with \code{\link[base]{cat}} for
#'   other repetitions of strings.
#' @family string repetitions
#' @examples
#'
#' hello(5) #prints the phrase 5 times
#'
#' @export
#'

hello <- function(times) {
  print(cat(rep("Hello, world!", times = times)))
}





