# Intro
## Project overview
Investigating role of collagen iii in models of virus-awakened dormant metastatic breast cancer.
## Directory
### QuPath
Project directory containing the main .`qpproj `file, `scripts`, image metadata and annotations in `data`, `classifiers`, and `resources`.

**Important:** image stacks not included, but filenames remain unaltered; when opening `project.qpproj`, point project at correct image directory and URIs should autofill
### RScripts
R analysis pipeline.
### data
Contains `.csv` files of raw cell counts.
### results
..
### Slides_Decoded_\[YYYYMMDD].csv
A .csv containing slide name and associated metadata. Included:

| Col | Feature        | Definition                                                                                                     |
| :-- | :------------- | -------------------------------------------------------------------------------------------------------------- |
| 1   | SlideName      | The slide ID as physically printed on the slide.                                                               |
| 2   | FileName       | The file as named. Typically matches the SlideName, with some exceptions.                                      |
| 3   | Model          | The mouse model from which the sample was derived                                                              |
| 4   | TrxGroup       | The virus treatment received by the model, if applicable                                                       |
| 5   | DPI            | Days post-infection, if applicable (see TrxGroup)                                                              |
| 6   | BatchStainDate | The date (YYYYMMDD) that the slide was stained.                                                                |
|     | BatchScanDate  | The date (YYYYMMDD) that the slide was scanned.                                                                |
| 7   | SampleNumber   | Refers to one of two samples on each slide. Either 1 (nearest to label) or 2 (further from label)              |
| 8   | Target1        | Target represented by channel 1: DAPI.                                                                         |
| 9   | Target2        | Target represented by channel 2: Her2                                                                          |
| 10  | Target3        | Target represented by channel 3: variable                                                                      |
| 11  | Channel1       | Fluorescence imaging channel 1: DAPI                                                                           |
| 12  | Channel2       | Fluorescence imaging channel 2: TRITC                                                                          |
| 13  | Channel3       | Fluorescence imaging channel 3: Cy5                                                                            |
| 14  | Exception      | Place for nonstandard comments, e.g. reason that the slide may be excluded from analysis, filename discrepancy |

# QuPath

> [!tip] DATE CONVENTION
> Dates always represented in YYYYMMDD in scripts, metadata, etc. for automation and consistency

## Scripts
### ...
# R Analysis

# Example Workflow
``` bash

```
