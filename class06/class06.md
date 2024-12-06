# Class 6: Functions
Emily Hendrickson (PID: A69034780)

**Functions!**

*fname \<- function(arg1, arg2) { paste(arg1,arg2)}*

My first function:

``` r
add <- function(x,y) {
  x + y
}

w <- c(1,3)

v <- c(2,4)

add(w,v)
```

    [1] 3 7

``` r
add(1,1)
```

    [1] 2

``` r
add(x=1, y=2)
```

    [1] 3

``` r
add(c(1,2,3), 3)
```

    [1] 4 5 6

To set a default value:

``` r
add <- function(x,y=5) {
  x + y
}

w <- c(1,3)

v <- c(2,4)

add(w,v)
```

    [1] 3 7

``` r
add(1,1)
```

    [1] 2

``` r
add(x=1, y=2)
```

    [1] 3

``` r
add(c(1,2,3), 3)
```

    [1] 4 5 6

``` r
add(10)
```

    [1] 15

**New Task:** To make a function *generatedna()* that makes a random
nucleotide sequence

``` r
generatedata <- function(length) {
  bases <- c("A","T","C","G")
  dna_sequence <- sample(x=bases, size = length, replace = TRUE)
  return(dna_sequence)
}

generatedata(5)
```

    [1] "C" "T" "A" "T" "T"

To make a function *generateprotein()* that makes a random amino acid
sequence.

``` r
generateprotein <- function(length) {
  residues <- unique(bio3d::aa.table$aa1)
  aa_sequence <- paste(sample(x=residues, size = length, replace = TRUE), 
                       collapse = "")
  return(aa_sequence)
}

# Generate a sequence of 5 amino acids
generateprotein(5)
```

    [1] "CMDKA"

To generate many sequences:

generate protein sequences of lengths 6 to 12 in fasta format

``` r
sequences <- sapply(6:12, generateprotein)

cat(paste(">id.length", 6:12, "\n", sequences, sep=""), sep="\n")
```

    >id.length6
    EWKPAE
    >id.length7
    YCCFAAQ
    >id.length8
    EYWCXDWA
    >id.length9
    SVAYSVQTF
    >id.length10
    RQHRNALMYH
    >id.length11
    CNSNCQPISKR
    >id.length12
    IXGXEEFGXEFD
