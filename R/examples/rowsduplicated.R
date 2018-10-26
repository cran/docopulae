a1 = as.double(seq(0, 0.4, 0.1))
a2 = as.double(c(0.0, 0.1, 0.2, 0.3, 0.4))

b1 = as.double(seq(-0.3, 0, 0.1))
b2 = as.double(c(-0.3, -0.2, -0.1, 0.0))

m1 = as.matrix(expand.grid(a1, b1))
m2 = as.matrix(expand.grid(a2, b2))

all(m1 == m2)  # probably FALSE due to precision of seq
all(as.character(m1) == as.character(m2))  # TRUE by design

m12 = rbind(m1, m2)

sum(duplicated(m12))  # depending on version of R, returns 20 or 8
if (sum(rowsduplicated(m12)) != 8)  # always returns 8
    stop('incorrect')
