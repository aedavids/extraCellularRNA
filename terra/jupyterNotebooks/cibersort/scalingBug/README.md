# cibersort scaling bug

When I implemented the DESeq2 normalization algorithm to test how well it worked with simulated mass spec data I discovered that we need to divide the samples by the estimated scaling factors not multiple

These notebooks implment division. I did not want to merge them back to main becuase the deconvolution fractions as k-way true positive results where worse than when we multipled by the scaling factors

My guess is the DESeq code inverts the scaling factors. TODO look at the source code to confirm my suspicion
