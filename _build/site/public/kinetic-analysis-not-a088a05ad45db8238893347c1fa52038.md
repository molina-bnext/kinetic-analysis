
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

:::{note}
Ideally, we would have a tool that (1) let's people describe experiments in a standardized way (spec lang), (2) that can then be mapped onto a platemap, and (3) displayed in part or whole in this section [this might be an interesting place for a widget to explore the parameters of the experiment]
:::

This is text describing the context of the experiment. There a few scenarios where this text can end up here. We can imagine users working in jupyterhub typing directly or its copy and pasted in from an experimental notebook. In this particular case, we have two experiments on testing the effect of DNA concentration and a second experiment comparing different DNA constrcuts - the infamous artifact 11 (AR-11)
 
:::{table} Experimental conditions studied
:label: pd:platemap
![](#pd:platemap)
:::

### Analysis

One thing that might be interesting to show here would be a progression from end-point analysis, unannotated timeseries, time series annotations, and then exploring the kinetic parameters. 

```{code} python
:label: code-block
:caption: we use this code to generate annotated figures
# select for a specific sub-experiment
data_concentration = data[data["Experiment"] == "DNA Concentration"]
```

The [](code-block) is used to generate the following figure using the [cytosol analysis toolkit](https://github.com/bnext-bio/nucleus/tree/main/cdk/analysis/cytosol-analysis)

:::{figure} #plt:expression-kinetics-simple
:label: fig:expression kinetics
:align: center
Correlation of top 1% liposome fluorescence to input DNA concentration. It appears that there may be two expression regimes, or a saturation effect.
:::

:::{figure} #plt:single-kinetics
:label: fig:kinetics-plot
:align: center
Correlation of top 1% liposome fluorescence to input DNA concentration. It appears that there may be two expression regimes, or a saturation effect.
:::
