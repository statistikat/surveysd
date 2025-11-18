# Methodology

In the following we present the methodology in `surveysd` by applying
the workflow described in
[`vignette("surveysd")`](https://statistikat.github.io/surveysd/articles/surveysd.md)
to multiple consecutive years of EU-SILC data for one country. The
methodology contains the following steps, in this order

- Draw $B$ bootstrap replicates from EU-SILC data for each year $y_{t}$,
  $t = 1,\ldots,n_{y}$ separately. Since EU-SILC has a rotating panel
  design the bootstrap replicate of a household is carried forward
  through the years. That is, the bootstrap replicate of a household in
  the follow-up years is set equal to the bootstrap replicate of the
  same household when it first enters EU-SILC.
- Multiply each set of bootstrap replicates by the sampling weights to
  obtain uncalibrated bootstrap weights and calibrate each of the
  uncalibrated bootstrap weights using iterative proportional fitting.
- Estimate the point estimate of interest $\theta$, for each year and
  each calibrated bootstrap weight to obtain
  ${\widetilde{\theta}}^{(i,y_{t})}$, $t = 1,\ldots,n_{y}$,
  $i = 1,\ldots,B$. For fixed $y_{t}$ apply a filter with equal weights
  for each $i$ on ${\widetilde{\theta}}^{(i,y^{*})}$,
  $y^{*} \in \{ y_{t - 1},y_{t},y_{t + 1}\}$ , to obtain
  ${\widetilde{\theta}}^{(i,y_{t})}$. Estimate the variance of $\theta$
  using the distribution of ${\widetilde{\theta}}^{(i,y_{t})}$.

## Bootstrapping

Bootstrapping has long been around and used widely to estimate
confidence intervals and standard errors of point estimates.\[Efron
(1979)} Given a random sample $\left( X_{1},\ldots,X_{n} \right)$ drawn
from an unknown distribution $F$ the distribution of a point estimate
$\theta\left( X_{1},\ldots,X_{n};F \right)$ can in many cases not be
determined analytically. However when using bootstrapping one can
simulate the distribution of $\theta$.

Let $s_{(.)}$ be a bootstrap sample, e.g. drawing $n$ observations with
replacement from the sample $\left( X_{1},\ldots,X_{n} \right)$, then
one can estimate the standard deviation of $\theta$ using $B$ bootstrap
samples through
$$sd(\theta) = \sqrt{\frac{1}{B - 1}\sum\limits_{i = 1}^{B}\left( \theta\left( s_{i} \right) - \overline{\theta} \right)^{2}}\quad,$$

with
$\overline{\theta}:=\frac{1}{B}\sum\limits_{i = 1}^{B}\theta\left( s_{i} \right)$
as the sample mean over all bootstrap samples.

In context of sample surveys with sampling weights one can use
bootstrapping to calculate so called bootstrap weights. These are
computed via the bootstrap samples $s_{i}$, $i = 1,\ldots,B$, where for
each $s_{i}$ every unit of the original sample can appear $0$- to
$m$-times. With $f_{j}^{i}$ as the frequency of occurrence of
observation $j$ in bootstrap sample $s_{i}$ the uncalibrated bootstrap
weights ${\widetilde{b}}_{j}^{i}$ are defined as:

$${\widetilde{b}}_{j}^{i} = f_{j}^{i}w_{j}\quad,$$

with $w_{j}$ as the calibrated sampling weight of the original sample.
Using iterative proportional fitting procedures one can recalibrate the
bootstrap weights ${\widetilde{b}}_{j}^{.}$, $j = 1,\ldots,B$ to get the
adapted or calibrated bootstrap weights $b_{j}^{i}$, $j = 1,\ldots,B$.

### Rescaled Bootstrap

Since EU-SILC is a stratified sample without replacement drawn from a
finite population the naive bootstrap procedure, as described above,
does not take into account the heterogeneous inclusion probabilities of
each sample unit. Thus it will not yield satisfactory results. Therefore
we will use the so called rescaled bootstrap procedure introduced and
investigated by (Rao and Wu 1988). The bootstrap samples are selected
without replacement and do incorporate the stratification as well as
clustering on multiple stages (see (Chipperfield and Preston
2007),(Preston 2009)).

For simplistic reasons we will only describe the rescaled bootstrap
procedure for a two stage stratified sampling design. For more details
on a general formulation please see (Preston 2009).

### Sampling design

Consider the finite population $U$ which is divided into $H$
non-overlapping strata $\bigcup\limits_{h = 1,\ldots,H}U_{h} = U$, of
which each strata $h$ contains of $N_{h}$ clusters. For each strata $h$,
$C_{hc}$, $c = 1,\ldots,n_{h}$ clusters are drawn, containing $N_{hc}$
households. Furthermore in each cluster $C_{hc}$ of each strata $h$
simple random sampling is performed to select a set of households
$Y_{hcj}$, $j = 1,\ldots,n_{hc}$.

### Bootstrap procedure

In contrast to the naive bootstrap procedure where for a stage,
containing $n$ sampling units, the bootstrap replicate is obtained by
drawing $n$ sampling units with replacement, for the rescaled bootstrap
procedure $n^{*} = \left\lfloor \frac{n}{2} \right\rfloor$ sampling
units are drawn without replacement. Given a value $x$,
$\lfloor x\rfloor$ denotes the largest integer smaller than $x$, whereas
$\lceil x\rceil$ denotes the smallest integer lager then $x$.
(Chipperfield and Preston 2007) have shown that the choice of either
$\left\lfloor \frac{n}{2} \right\rfloor$ or
$\left\lceil \frac{n}{2} \right\rceil$ is optimal for bootstrap samples
without replacement, although $\left\lfloor \frac{n}{2} \right\rfloor$
has the desirable property that the resulting uncalibrated bootstrap
weights will never be negative.

At the first stage the $i$-th bootstrap replicate, $f_{hc}^{i,1}$, for
each cluster $C_{hc}$,$c = 1,\ldots,n_{h}$, belonging to strata $h$, is
defined by

$$f_{hc}^{i,1} = 1 - \lambda_{h} + \lambda_{h}\frac{n_{h}}{n_{h}^{*}}\delta_{hc}\quad\quad\forall c \in \{ 1,\ldots,n_{h}\}$$
with
$$n_{h}^{*} = \left\lfloor \frac{n_{h}}{2} \right\rfloor$$$$\lambda_{h} = \sqrt{\frac{n_{h}^{*}\left( 1 - \frac{n_{h}}{N_{h}} \right)}{n_{h} - n_{h}^{*}}}\quad,$$

where $\delta_{hc} = 1$ if cluster $c$ is selected in the sub-sample of
size $n_{h}^{*}$ and 0 otherwise.

The $i$-th bootstrap replicate at the second stage, $f_{hcj}^{i,2}$, for
each household $Y_{hcj}$, $j = 1,\ldots,n_{hc}$, belonging to cluster
$c$ in strata $h$ is defined by

$$f_{hcj}^{i,2} = f_{hc}^{i,1} - \lambda_{hc}\sqrt{\frac{n_{h}}{n_{h}^{*}}}\delta_{hc}\left\lbrack \frac{n_{hc}}{n_{hc}^{*}}\delta_{hcj} - 1 \right\rbrack\quad\quad\forall c \in \{ 1,\ldots,n_{h}\}$$
with
$$n_{hc}^{*} = \left\lfloor \frac{n_{hc}}{2} \right\rfloor$$$$\lambda_{hc} = \sqrt{\frac{n_{hc}^{*}N_{h}\left( 1 - \frac{n_{hc}}{N_{hc}} \right)}{n_{hc} - n_{hc}^{*}}}\quad,$$

where $\delta_{hcj} = 1$ if household $j$ is selected in the sub sample
of size $n_{hc}^{*}$ and 0 otherwise.

### Single PSUs

When dealing with multistage sampling designs the issue of single PSUs,
e.g. a single response unit is present at a stage or in a strata, can
occur. When applying bootstrapping procedures these single PSUs can lead
to a variety of issues. For the methodology proposed in this work we
combined single PSUs at each stage with the next smallest strata or
cluster, before applying the bootstrap procedure.

### Taking bootstrap replicates forward

The bootstrap procedure above is applied on the EU-SILC data for each
year $y_{t}$, $t = 1,\ldots,n_{y}$ separately. Since EU-SILC is a yearly
survey with rotating penal design the $i$-th bootstrap replicate at the
second stage, $f_{hcj}^{i,2}$, for a household $Y_{hcj}$ is taken
forward until the household $Y_{hcj}$ drops out of the sample. That is,
for the household $Y_{hcj}$, which enters EU-SILC at year $y_{1}$ and
drops out at year $y_{\widetilde{t}}$, the bootstrap replicates for the
years $y_{2},\ldots,y_{\widetilde{t}}$ are set to the bootstrap
replicate of the year $y_{1}$.

### Split households

Due to the rotating penal design so called split households can occur.
For a household participating in the EU-SILC survey it is possible that
one or more residents move to a new so called split household, which is
followed up on in the next wave. To take this dynamic into account we
extended the procedure of taking forward the bootstrap replicate of a
household for consecutive waves of EU-SILC by taking forward the
bootstrap replicate to the split household. That means, that also any
new individuals in the split household will inherit this bootstrap
replicate.

Taking bootstrap replicates forward as well as considering split
households ensures that bootstrap replicates are more comparable in
structure with the actual design of EU-SILC.

### Uncalibrated bootstrap weights

Using the $i$-th bootstrap replicates at the second stage one can
calculate the $i$-th uncalibrated bootstrap weights $b_{hcj}^{i}$ for
each household $Y_{hcj}$ in cluster $c$ contained in strata $h$ by

$${\widetilde{b}}_{hcj}^{i} = f_{hcj}^{i,2}w_{hcj}\quad,$$ where
$w_{hcj}$ corresponds to the original household weight contained in the
sample.

For ease of readability we will drop the subindices regarding strata $h$
and cluster $c$ for the following sections, meaning that the $j$-th
household in cluster $c$ contained in strata $h$, $Y_{hcj}$, will now be
denoted as the $j$-th household, $Y_{j}$, where $j$ is the position of
the household in the data. In accordance to this the $i$-th uncalibrated
bootstrap replicates for household $j$ are thus denoted as
${\widetilde{b}}_{j}^{i}$ and the original household weight as $w_{j}$.

## Iterative proportional fitting (IPF)

The uncalibrated bootstrap weights ${\widetilde{b}}_{j}^{i}$ computed
through the rescaled bootstrap procedure yields population statistics
that differ from the known population margins of specified
sociodemographic variables for which the base weights $w_{j}$ have been
calibrated. To adjust for this the bootstrap weights
${\widetilde{b}}_{j}^{i}$ can be recalibrated using iterative
proportional fitting as described in (Meraner, Gumprecht, and Kowarik
2016).

Let the original weight $w_{j}$ be calibrated for $n = n_{P} + n_{H}$
sociodemographic variables which are divided into the sets
$\mathcal{P}:=\{ p_{c},c = 1\ldots,n_{P}\}$ and
$\mathcal{H}:=\{ h_{c},c = 1\ldots,n_{H}\}$. $\mathcal{P}$ and
$\mathcal{H}$ correspond to personal, for example gender or age, or
household variables, like region or households size, respectively. Each
variable in either $\mathcal{P}$ or $\mathcal{H}$ can take on $P_{c}$ or
$H_{c}$ values with and $N_{v}^{p_{c}}$, $v = 1,\ldots,P_{c}$, or
$N_{v}^{h_{c}}$, $v = 1,\ldots,H_{c}$, as the corresponding population
margins. Starting with $k = 0$ the iterative proportional fitting
procedure is applied on each ${\widetilde{b}}_{j}^{i}$, $i = 1,\ldots,B$
separately. The weights are first updated for personal and afterwards
updated for household variables. If constraints regarding the
populations margins are not met $k$ is raised by 1 and the procedure
starts from the beginning. For the following denote as starting weight
${\widetilde{b}}_{j}^{\lbrack 0\rbrack}:={\widetilde{b}}_{j}^{i}$ for
fixed $i$.

### Adjustment and trimming for $\mathcal{P}$

The uncalibrated bootstrap weight
${\widetilde{b}}_{j}^{\lbrack{(n + 1)}k + c - 1\rbrack}$ for the $j$-th
observation is iteratively multiplied by a factor so that the projected
distribution of the population matches the respective calibration
specification $N_{p_{c}}$, $c = 1,\ldots,n_{P}$. For each
$c \in \left\{ 1,\ldots,n_{P} \right\}$ the calibrated weights against
$N_{v}^{p_{c}}$ are computed as
$${\widetilde{b}}_{j}^{\lbrack{(n + 1)}k + c\rbrack} = {{\widetilde{b}}_{j}}^{\lbrack{(n + 1)}k + c - 1\rbrack}\frac{N_{v}^{p_{c}}}{\sum\limits_{l}{\widetilde{b}}_{l}^{\lbrack{(n + 1)}k + c - 1\rbrack}},$$
where the summation in the denominator expands over all observations
which have the same value as observation $j$ for the sociodemographic
variable $p_{c}$. If any weights
${\widetilde{b}}_{j}^{\lbrack nk + c\rbrack}$ fall outside the range
$\left\lbrack \frac{w_{j}}{4};4w_{j} \right\rbrack$ they will be recoded
to the nearest of the two boundaries. The choice of the boundaries
results from expert-based opinions and restricts the variance of which
has a positive effect on the sampling error. This procedure represents a
common form of weight trimming where very large or small weights are
trimmed in order to reduce variance in exchange for a possible increase
in bias ((Potter 1990),(Potter 1993)).

### Averaging weights within households

Since the sociodemographic variables $p_{1},\ldots,p_{n_{c}}$ include
person-specific variables, the weights
${\widetilde{b}}_{j}^{\lbrack nk + n_{p}\rbrack}$ resulting from the
iterative multiplication can be unequal for members of the same
household. This can lead to inconsistencies between results projected
with household and person weights. To avoid such inconsistencies each
household member is assigned the mean of the household weights. That is
for each person $j$ in household $a$ with $h_{a}$ household members, the
weights are defined by
$${\widetilde{b}}_{j}^{\lbrack{(n + 1)}k + n_{p} + 1\rbrack} = \frac{\sum\limits_{l \in a}{\widetilde{b}}_{l}^{\lbrack{(n + 1)}k + n_{p}\rbrack}}{h_{a}}$$
This can result in losing the population structure performed in the
previous subsection.

### Adjustment and trimming for $\mathcal{H}$

After adjustment for individual variables the weights
$b_{j}^{\lbrack nk + n_{p} + 1\rbrack}$ are updated for the set of
household variables $\mathcal{H}$ according to a household convergence
constraint parameter $\epsilon_{h}$. The parameters $\epsilon_{h}$
represent the allowed deviation from the population margins using the
weights $b_{j}^{\lbrack nk + n_{p} + 1\rbrack}$ compared to
$N_{v}^{h_{c}}$, $c = 1,\ldots,n_{H}$, $v = 1,\ldots,H_{c}$. The updated
weights are computed as
$$b_{j}^{\lbrack{(n + 1)}k + n_{p} + c + 1\rbrack} = \left\{ \begin{array}{l}
{b_{j}^{\lbrack{(n + 1)}k + n_{p} + 1\rbrack}\frac{N_{v}^{h_{c}}}{\sum\limits_{l}b_{l}^{\lbrack{(n + 1)}k + n_{p} + 1\rbrack}}\quad{\text{if}\mspace{6mu}}\sum\limits_{l}b_{j}^{\lbrack{(n + 1)}k + n_{p} + 1\rbrack} \notin \left( \left( 1 - 0.9\epsilon_{h} \right)N_{v}^{h_{c}},\left( 1 + 0.9\epsilon_{h} \right)N_{v}^{h_{c}} \right)} \\
{b_{j}^{\lbrack{(n + 1)}k + n_{p} + 1\rbrack}\quad\text{otherwise}}
\end{array} \right.$$ with the summation in the denominator ranging over
all households $l$ which take on the same values for $h_{c}$ as
observation $j$. As described in the previous subsection the new weight
are recoded if they exceed the interval
$\left\lbrack \frac{w_{j}}{4};4w_{j} \right\rbrack$ and set to the upper
or lower bound, depending of
$b_{j}^{\lbrack{(n + 1)}k + n_{p} + c + 1\rbrack}$ falls below or above
the interval respectively.

### Convergence

For each adjustment and trimming step the factor
$\frac{N_{v}^{(.)}}{\sum\limits_{l}b_{l}^{\lbrack{(n + 1)}k + j\rbrack}}$,
$j \in \{ 1,\ldots,n + 1\} \smallsetminus \{ n_{p} + 1\}$, is checked
against convergence constraints for households, $\epsilon_{h}$, or
personal variables $\epsilon_{p}$, where $(.)$ corresponds to either a
household or personal variable. To be more precise for variables in
$\mathcal{P}$ the constraints

$$\frac{N_{v}^{p_{c}}}{\sum\limits_{l}{\widetilde{b}}_{l}^{\lbrack{(n + 1)}k + j\rbrack}} \in \left( \left( 1 - \epsilon_{p} \right)N_{v}^{p_{c}},\left( 1 + \epsilon_{p} \right)N_{v}^{p_{c}} \right)$$
and for variables in $\mathcal{H}$ the constraints

$$\frac{N_{v}^{h_{c}}}{\sum\limits_{l}{\widetilde{b}}_{l}^{\lbrack{(n + 1)}k + j\rbrack}} \in \left( \left( 1 - \epsilon_{h} \right)N_{v}^{h_{c}},\left( 1 + \epsilon_{h} \right)N_{v}^{h_{c}} \right)$$
are verified, where the sum in the denominator expands over all
observations which have the same value for variables $h_{c}$ or $p_{c}$.
If these constraints hold true the algorithm reaches convergence,
otherwise $k$ is raised by 1 and the procedure repeats itself.

The above described calibration procedure is applied on each year
$y_{t}$ of EU-SILC separately, $t = 1,\ldots n_{y}$, thus resulting in
so called calibrated bootstrap sample weights $b_{j}^{(i,y_{t})}$,
$i = 1,\ldots,B$ for each year $y$ and each household $j$.

## Variance estimation

Applying the previously described algorithms to EU-SILC data for
multiple consecutive years $y_{t}$, $t = 1,\ldots n_{y}$, yields
calibrated bootstrap sample weights $b_{j}^{(i,y_{t})}$ for each year
$y_{t}$. Using the calibrated bootstrap sample weights it is straight
forward to compute the standard error of a point estimate
$\theta\left( \textbf{𝐗}^{y_{t}},\textbf{𝐰}^{y_{t}} \right)$ for year
$y_{t}$ with
$\textbf{𝐗}^{y_{t}} = \left( X_{1}^{y_{t}},\ldots,X_{n}^{y_{t}} \right)$
as the vector of observations for the variable of interest in the survey
and $\textbf{𝐰}^{y_{t}} = (w_{1}^{y_{t}},\ldots,w_{n}^{y_{t}}$ as the
corresponding weight vector, with

$$sd(\theta) = \sqrt{\frac{1}{B - 1}\sum\limits_{i = 1}^{B}\left( \theta^{(i,y_{t})} - \overline{\theta^{(.,y_{t})}} \right)^{2}}$$
with
$$\overline{\theta^{(.,y_{t})}} = \frac{1}{B}\sum\limits_{i = 1}^{B}\theta^{(i,y_{t})}\quad,$$
where
$\theta^{(i,y_{t})}:=\theta\left( \textbf{𝐗}^{y_{t}},\textbf{𝐛}^{(i,y_{t})} \right)$
is the estimate of $\theta$ in the year $y_{t}$ using the $i$-th vector
of calibrated bootstrap weights.

As already mentioned the standard error estimation for indicators in
EU-SILC yields high quality results for NUTS1 or country level. When
estimation indicators on regional or other sub-aggregate levels one is
confronted with point estimates yielding high variance.

To overcome this issue we propose to estimate $\theta$ for 3,
consecutive years using the calibrated bootstrap weights, thus
calculating
$\{\theta^{(i,y_{t - 1})},\theta^{(i,y_{t})},\theta^{(i,y_{t + 1})}\}$,
$i = 1,\ldots,B$. For fixed $i$ one can apply a filter with equal filter
weights on the time series
$\{\theta^{(i,y_{t - 1})},\theta^{(i,y_{t})},\theta^{(i,y_{t + 1})}\}$
to create ${\widetilde{\theta}}^{(i,y_{t})}$

$${\widetilde{\theta}}^{(i,y_{t})} = \frac{1}{3}\left\lbrack \theta^{(i,y_{t - 1})} + \theta^{(i,y_{t})} + \theta^{(i,y_{t + 1})} \right\rbrack\quad.$$

Doing this for all $i$, $i = 1,\ldots,B$, yields
${\widetilde{\theta}}^{(i,y_{t})}$, $i = 1,\ldots,B$. The standard error
of $\theta$ can then be estimated with

$$sd(\theta) = \sqrt{\frac{1}{B - 1}\sum\limits_{i = 1}^{B}\left( {\widetilde{\theta}}^{(i,y_{t})} - \overline{{\widetilde{\theta}}^{(.,y_{t})}} \right)^{2}}$$
with
$$\overline{{\widetilde{\theta}}^{(.,y_{t})}} = \frac{1}{B}\sum\limits_{i = 1}^{B}{\widetilde{\theta}}^{(i,y_{t})}\quad.$$

Applying the filter over the time series of estimated
$\theta^{(i,y_{t})}$ leads to a reduction of variance for $\theta$ since
the filter reduces the noise in
$\{\theta^{(i,y_{t - 1})},\theta^{(i,y_{t})},\theta^{(i,y_{t + 1})}\}$
and thus leading to a more narrow distribution for
${\widetilde{\theta}}^{(i,y_{t})}$.

It should also be noted that estimating indicators from a survey with
rotating panel design is in general not straight forward because of the
high correlation between consecutive years. However with our approach to
use bootstrap weights, which are independent from each other, we can
bypass the cumbersome calculation of various correlations, and apply
them directly to estimate the standard error. (Bauer et al. 2013) showed
that using the proposed method on EU-SILC data for Austria the reduction
in resulting standard errors corresponds in a theoretical increase in
sample size by about 25$\%$. Furthermore this study compared this method
to the use of small area estimation techniques and on average the use of
bootstrap sample weights yielded more stable results.

## References

Bauer, Martin, Matthias Till, Richard Heuberger, Marcel Bilgili, Thomas
Glaser, Elisabeth Kafka, Johannes Klotz, et al. 2013. “Studie Zu Armut
Und Sozialer Eingliederung in Den Bundesl"andern.” Statistik Austria
\[in German\].

Chipperfield, James, and John Preston. 2007. “Efficient Bootstrap for
Business Surveys.” *Survey Methodology* 33 (December): 167–72.
<https://www150.statcan.gc.ca/n1/en/catalogue/12-001-X200700210494>.

Efron, B. 1979. “Bootstrap Methods: Another Look at the Jackknife.”
*Ann. Statist.* 7 (1): 1–26. <https://doi.org/10.1214/aos/1176344552>.

Meraner, Angelika, Daniela Gumprecht, and Alexander Kowarik. 2016.
“Weighting Procedure of the Austrian Microcensus Using Administrative
Data.” *Austrian Journal of Statistics* 45 (June): 3.
<https://doi.org/10.17713/ajs.v45i3.120>.

Potter, Frank J. 1990. “A Study of Procedures to Identify and Trim
Extreme Sampling Weights.” *Proceedings of the American Statistical
Association, Section on Survey Research Methods*, 225–30.
<http://www.asasrms.org/Proceedings/papers/1990_034.pdf>.

———. 1993. “The Effect of Weight Trimming on Nonlinear Survey
Estimates.” *Proceedings of the American Statistical Association,
Section on Survey Research Methods* 2: 758–63.
<http://www.asasrms.org/Proceedings/papers/1993_127.pdf>.

Preston, J. 2009. “Rescaled Bootstrap for Stratified Multistage
Sampling.” *Survey Methodology* 35 (December): 227–34.
<https://www150.statcan.gc.ca/n1/en/catalogue/12-001-X200900211044>.

Rao, J. N. K., and C. F. J. Wu. 1988. “Resampling Inference with Complex
Survey Data.” *Journal of the American Statistical Association* 83
(401): 231–41.
