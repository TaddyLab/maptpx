map estimation of topic models
======

This is the maptpx package for R.  It implements map estimation for latent topic models in text analysis, as described in my paper <a href="http://jmlr.csail.mit.edu/proceedings/papers/v22/taddy12/taddy12.pdf">"On Estimation and Selection for Topic Models"</a>.  These programs were previously part of <a href="http://www.cran.r-project.org/web/packages/textir/index.html">textir</a>, but were spun off at version 2.0 of that package.  maptpx is no longer actively maintained; instead, we are focusing on fast distributed factor models based on <a href="http://arxiv.org/abs/1311.6139">distributed multinomial regression</a>, as implemented in the <a href="http://cran.r-project.org/web/packages/distrom/index.html">distrom</a> package.

If you want to take advantage of openmp parallelization, uncomment the relevant flags in src/MAKEVARS before compiling.
