# Rao-Wu Bootstrap in surveysd

In the following we describe an extension of the methodology implemented
in `surveysd`. While the default rescaling method described in
`vignette("Methodology")` follows the approach of Preston, this vignette
introduces an alternative rescaling method based on the Rao-Wu bootstrap
(Rao and Wu 1988).

### Mathematical Formulation

More formally, the replicate weights are constructed as follows:
$$w_{hi}^{*} = \left( 1 - \lambda_{h} + \lambda_{h} \cdot \frac{n_{h}}{m_{h}} \cdot r_{hi}^{*} \right) \cdot w_{hi}$$
with
$$\lambda_{h} = \sqrt{\frac{m_{h}\left( 1 - f_{h} \right)}{n_{h} - 1}}$$
where:  
- $w_{hi}$ is the design weight for PSU $i$ in stratum $h$,
$w_{h} = N_{h}/n_{h}$ is the average design weight for the entire
stratum $h$,  
- $N_{h}$ is the total number of units in stratum $h$,  
- $n_{h}$ is the number of PSUs in the original sample for stratum
$h$,  
- $r_{hi}^{*} \in \{ 0,1,2,\ldots\}$ is a **resampling indicator** that
represents how many times PSU $i$ in stratum $h$ was drawn in the
bootstrap ,  
- $m_{h} = n_{h} - 1$ is the **number of units** to be drawn in each
replicate for stratum $h$,  
- $\lambda_{h}$ is the **scaling factor** used to adjust the weights
during the bootstrap process,  
- and $f_{h} = n_{h}/N_{h}$ is the **sampling fraction** in stratum
$h$.  

### Use of Rao-Wu

The Rao-Wu bootstrap is particularly suited for complex survey designs
with stratification and multi-stage selection. The design-based variance
estimation benefits from a resampling method that accounts for the
primary sampling stage and incorporates finite population corrections
(FPC) directly in the replicate weights. The method is appropriate when
the sampling fraction in the first stage (i.e., the share of selected
PSUs within strata) is relatively small, typically below 10%.

Rao-Wu is **not suitable** when:  
- The sampling fraction in the first stage is large and the first-stage
sampling is without replacement.  
- The survey does not include PSU-level identifiers, making a
first-stage resampling infeasible since the method requires resampling
entire PSUs to reflect first-stage sampling variability.  
- The design is single-stage, where methods like the Preston bootstrap
may be more appropriate.  

In those cases, it is recommended to fall back on alternative bootstrap
methods such as the default Preston approach implemented in `surveysd`,
which offers more flexibility for a wider range of designs without PSU
information.

### Single PSUs

When dealing with multistage sampling designs, the issue of single PSUs,
e.g. a single response unit at a stage or in a stratum, can arise. When
applying resampling methods such as the Rao-Wu bootstrap, these single
PSUs can introduce challenges in the resampling process. In the Rao-Wu
method, we adjust for the presence of single PSUs by combining them with
the next smallest stratum or cluster before applying the resampling
procedure. This ensures that the resampling reflects the structure of
the original design while maintaining appropriate variance estimates for
the total and other statistics of interest.

### Example Implementation

#### Load Dataset

``` r
library(surveysd)

set.seed(1234)
eusilc <- demo.eusilc(n = 2, prettyNames = TRUE)

eusilc[1:5, .(year, povertyRisk, gender, pWeight)]
```

#### Draw bootstrap replicates

For the bootstrap select ‘method = “Rao-Wu”’. Otherwise the default
“Preston” is used.

``` r
dat_boot_rw <- draw.bootstrap(eusilc, 
                              method = "Rao-Wu",
                              REP = 10, 
                              hid = "hid", 
                              weights = "pWeight", 
                              strata = "region", 
                              period = "year")
```

#### Calibrate bootstrap replicates

Calibrate each sample according to the distribution of `gender` (on a
personal level) and `region` (on a household level).

``` r
dat_boot_calib <- recalib(dat_boot_rw, 
                          conP.var = "gender", 
                          conH.var = "region",
                          epsP = 1e-2, 
                          epsH = 2.5e-2, 
                          verbose = FALSE)
dat_boot_calib[1:5, .(year, povertyRisk, gender, pWeight, w1, w2, w3, w4)]
```

Rao, J. N. K., and C. F. J. Wu. 1988. “Resampling Inference with Complex
Survey Data.” *Journal of the American Statistical Association* 83
(401): 231–41.
