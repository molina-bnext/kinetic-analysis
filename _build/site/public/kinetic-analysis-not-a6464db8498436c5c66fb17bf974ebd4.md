
# Time course measurements

A shared understanding of how to interpret data from CFE systems would speed the development of measurements and standards towards improved reproducibility and address challenges in interpreting and comparing existing and future data, including data from failed experiments. Workshop participants agreed that time-course measurements of the product expressed in the CFE reaction should be favored over endpoint measurements, despite the time, labor, and costs involved, to obtain a more complete and informative view of a CFE reaction. Whenever possible, measurements should be reported as a reduced quantity, such as a mean value, with uncertainty (Figure 2A) and include a baseline from negative control measurements. Several recent studies have commented on metrics available from a time-course measurement of protein expression in a CFE workflow, for example, for the purpose of reaction optimization [ref] and the development of predictive modeling tools [ref]. Workshop participants emphasized several of these metrics, displayed below in Figure 2


- Maximum yield of product expressed;
- $t_{max}$ Time to reach the maximum yield of product expressed, as the time from the start of the measurement to the time to reach the maximum yield;
- $ v_{max}$ Maximum rate of product expression, as the maximum linear rate of production;
- Lag time, as the time from the start of the measurement to the time to reach the maximum rate of expression;
- Inflection time, as the time from the start of the measurement to the time at which the rate of product expression begins to decrease;
- Percent decline, as a decrease in the amount of product expressed from the maximum yield; and,
- Time to percent decline, as the amount of time from the start of the measurement to the time to reach a predefined decrease in the yield of product expressed after having reached its maximum yield of product expressed.

## Workflow

### Platemap


 
:::{table} Experimental conditions studied
:label: pd:platemap
![](#pd:platemap)
:::

### Analysis

One thing that might be interesting to show here would be a progression from end-point analysis, unannotated timeseries, time series annotations, and then exploring the kinetic parameters. 

```{code} python
:label: code:filter-platemap
:caption: filter the plate map to select for a specific experiment
# experiment -> "DNA Concentration"
data_concentration = data[data["Experiment"] == "DNA Concentration"]
```

Here, [](#code:filter-platemap) is used to generate the following figure using these [specific lines](https://github.com/bnext-bio/nucleus/blob/main/cdk/analysis/cytosol-analysis/cytosol-kinetics.ipynb#L563-L565) of code from the cytosol analysis toolkit. This might be useful for instance if we're showcasing a new feature in the codebase and we want to make a tight connection between code and function, inviting engagement.

The immediate advantage of this notebook over our existing way of publishing is that the figure is unambiguously linked to the data and code used to generate it. In this particular case, the sample data is co-located in the directory used to generate the HTML. Another approach might be to pull the data programmatically from an archive like zenodo.

:::{admonition} Pulling data directly from an archive
:class: dropdown
```{code} python
# this workflow would also work for OSF
from zenodo_client import Zenodo

zenodo = Zenodo(access_token=token)
OOH_NA_NA_RECORD = '13852103'

platemap_file = '20240916-platemap.xlsx'
data_file = 'timecourse.csv'

data = zenodo.download_latest(OOH_NA_NA_RECORD, data_file)
platemap = zenodo.download_latest(OOH_NA_NA_RECORD, platemap_file)
```
:::

:::{figure} #plt:expression-kinetics-simple
:label: fig:expression kinetics
:align: center
The plasmid pT7-deGFP (AR-11) is swept from 0-100 ng/uL in the PURE system
:::



:::{figure} #plt:single-kinetics
:label: fig:kinetics-plot
:align: center
Correlation of top 1% liposome fluorescence to input DNA concentration. It appears that there may be two expression regimes, or a saturation effect.
:::

