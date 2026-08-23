# Nice Printing of Multiline Texts

Useful for `print(fit$jags_code)`, `print(mcp_demo$example_code)`, etc.

## Usage

``` r
# S3 method for class 'mcptext'
print(x, ...)
```

## Arguments

- x:

  Character, often with newlines.

- ...:

  Currently ignored.

## Author

Jonas Kristoffer Lindeløv <jonas@lindeloev.dk>

## Examples

``` r
mytext = "line1 = 2\n line2 = 'horse'"
class(mytext) = "mcptext"
print(mytext)
#> line1 = 2
#>  line2 = 'horse'
```
