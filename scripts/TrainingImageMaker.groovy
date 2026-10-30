///////////////////////

// this script:
// 1. outlines your sample using defined tissueClassifier
// 2. creates a specified number of Region* annotations to create a training image
// 3. SKIPS any entries that do not match specified metadata
// run this then RemoveRegionAnnos to clean out Region* annos
// by Bailey 8/7/2026

//////////////////////

// set up to run over project each time
//def tissueClassifier =  "" // pixel classifier to select sample, leave commented out if you don't want to manually specify
def specifiedKey = "BatchScanDate" // what key you want to sort on (batch)
def specifiedValue = "20260424" // what value you want to match to
def tissueClassifier = "TisSlct_${specifiedKey}_${specifiedValue}"  //uncomment if you want to automatically set given Key and Value
def num_ROIs = 2 // how many ROIs per image
// dims (microns) of ROIs
double ROI_x_um = 500
double ROI_y_um = 500


// dont touch below this line

// Look for a specific metadata key
def entry = getProjectEntry()
def entryValue = entry.getMetadataValue(specifiedKey)

// Skip if the metadata doesn't match target
if ((entryValue != specifiedValue)||(entry.getMetadataValue("Type") != "Sample")||(entry.getMetadataValue("Cy5Target") != "Col3a1")) { // also excludes anything that's not a sample; in this case, samples have metadata type=tissue
    print "Skipping ${entry.getImageName()}"
    return 
}

print "Creating region annos for image ${entry.getImageName()}"


// set ROI size given avg pixel size
double micronsPerPixel = getCurrentImageData().getServer().getPixelCalibration().getAveragedPixelSize()
double ROI_x_px = ROI_x_um/micronsPerPixel // convert width to px
double ROI_y_px = ROI_y_um/micronsPerPixel // ^^ height


// select valid sample area
resetSelection()
createAnnotationsFromPixelClassifier(tissueClassifier, 100000.0, 5000000.0) // run tissue selector

// next: random rectangles script from Krstin Gallik
//setup some classes that will be used later and selecting the main annotation that the random rectangles will be in
def annoList = getAnnotationObjects().findAll{it.getPathClass() == getPathClass("Tissue")}
def anno = annoList.get(0)
def tissue = anno.getROI()
def tissueClass = getPathClass("Tissue")
def roiClass = getPathClass("Region*")

anno.each{it.setPathClass(tissueClass)} //assigned a classification to the main annotation so it is not included in future steps that are not shown in this code snippet
def random_ROIs = [] //empty list to store in the random rectangles

//basic for loop that creates the number of random rectangles calculated above
for (i in 1..num_ROIs) {
    randomROI = RoiTools.createRandomRectangle(tissue, ROI_x_px, ROI_y_px)
    randomAnnot = PathObjects.createAnnotationObject(randomROI)
    random_ROIs.add(randomAnnot)
    tissue = RoiTools.difference(tissue, randomROI) // this step is essentially punching holes into the original ROI and then is updating the tissue annotation being used to place the next ROI at the top of the loop. It keeps the newly added ROIs from overlapping eachother.
    //println tissue.getArea() //used this line when testing to verify that the main annotation was actually decreasing in size
}
// Add in the new ROIs and give them a class for easier selection should you want to do something with them
addObjects(random_ROIs)
def ROIs = getAnnotationObjects().findAll{it.getPathClass() == null}
ROIs.each{it.setPathClass(roiClass)}
// Once the ROIs are created, they can be used to grab measurements, or create detections within (what I used mine for)
// end Gallik script

