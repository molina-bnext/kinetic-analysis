
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

the general workflow will involves importing data of an appropriate format. 

| Time            |   T° | Well   |   Data | Row   |   Column | Read    |   Ex |   Em |   Seconds | Experiment        | Label              | DNA Template                |   [DNA Template] (ng/uL) | Blank   |   DataBlank |   DataBlanked |
|:----------------|-----:|:-------|-------:|:------|---------:|:--------|-----:|-----:|----------:|:------------------|:-------------------|:----------------------------|-------------------------:|:--------|------------:|--------------:|
| 0 days 00:00:00 | 37   | B1     |    289 | B     |        1 | 490,520 |  490 |  520 |         0 | DNA Concentration | GFP-100            | pT7-deGFP (AR-11)           |                    100   | nan     |         258 |            31 |
| 0 days 00:05:00 | 37.1 | B1     |    312 | B     |        1 | 490,520 |  490 |  520 |       300 | DNA Concentration | GFP-100            | pT7-deGFP (AR-11)           |                    100   | nan     |         247 |            65 |
| 0 days 00:10:00 | 37   | B1     |    449 | B     |        1 | 490,520 |  490 |  520 |       600 | DNA Concentration | GFP-100            | pT7-deGFP (AR-11)           |                    100   | nan     |         253 |           196 |
| 0 days 00:15:00 | 37   | B1     |    982 | B     |        1 | 490,520 |  490 |  520 |       900 | DNA Concentration | GFP-100            | pT7-deGFP (AR-11)           |                    100   | nan     |         246 |           736 |
| 0 days 00:20:00 | 37   | B1     |   1880 | B     |        1 | 490,520 |  490 |  520 |      1200 | DNA Concentration | GFP-100            | pT7-deGFP (AR-11)           |                    100   | nan     |         243 |          1637 |
| 0 days 00:25:00 | 37   | B1     |   3161 | B     |        1 | 490,520 |  490 |  520 |      1500 | DNA Concentration | GFP-100            | pT7-deGFP (AR-11)           |                    100   | nan     |         250 |          2911 |
| 0 days 00:30:00 | 37   | B1     |   4846 | B     |        1 | 490,520 |  490 |  520 |      1800 | DNA Concentration | GFP-100            | pT7-deGFP (AR-11)           |                    100   | nan     |         240 |          4606 |
| 0 days 00:35:00 | 36.9 | B1     |   6747 | B     |        1 | 490,520 |  490 |  520 |      2100 | DNA Concentration | GFP-100            | pT7-deGFP (AR-11)           |                    100   | nan     |         241 |          6506 |
| 0 days 00:40:00 | 37   | B1     |   8842 | B     |        1 | 490,520 |  490 |  520 |      2400 | DNA Concentration | GFP-100            | pT7-deGFP (AR-11)           |                    100   | nan     |         241 |          8601 |

```{code} python
:label: code-block
:caption: we use this code to generate annotated figures
import numpy as np

np.sin(x)
```

The [](code-block) is used to generate the following figure using the [cytosol analysis toolkit](https://github.com/bnext-bio/nucleus/tree/main/cdk/analysis/cytosol-analysis)

:::{figure} #plt:expression-panel-collapsed
:label: fig:mean-gfp-intensity-correlation
:align: center
Correlation of top 1% liposome fluorescence to input DNA concentration. It appears that there may be two expression regimes, or a saturation effect.
:::
