///////////////////////

// this script:
// 1. outlines your sample using defined tissueClassifier (if such an obj doesn't already exist)
// 2. removes any regions with specific ignoreClass
// 3. creates random ROIs over specified coverage
// 4. runs cell detector in random ROIs
// 5. runs specified cellClassifier
// 6. repeats if u want
// 7. SKIPS any entries that do not match specified metadata
// by Bailey 8/7/2026
// see end for troubleshooting

//////////////////////

// first: setup params
def entry = getProjectEntry() // leave this alone
def specifiedKey = "BatchScanDate" // what key you want to sort on (batch)

// CHOOSE ONE OR THE OTHER: AUTOMATIC BATCH PROCESSING OR MANUAL SINGLE-BATCH PROCESSING
    def specifiedValue = "20260420" // MANUALLY SET. what value you want to match to
    // def specifiedValue = entry.getMetadataValue(specifiedKey) // AUTOMATICALLY SET. uncomment if you want to match tissue selector and cell classifier to key value
    
// CHOOSE ONE OR THE OTHER: tissue classifier 
    def tissueClassifier = "TisSlct_${specifiedKey}_${specifiedValue}"  // AUTOMATICALLY SET. uncomment if you want to automatically set given Key and Value
    // may have to change prefix ^ here
    //def tissueClassifier = "SampleIDer 1.0" // MANUALLY SET the pixel classifier used to outline your sample

// CHOOSE ONE OR THE OTHER: cell classifier
    def cellClassifier = "Her2Dtct_${specifiedKey}_${specifiedValue}" // AUTOMATICALLY SET given metadata, uncomment if you want to automatically set 
    // may have to change prefix ^ here
    //def cellClassifier = "Her2+ 1.0" // MANUALLY SET name of cell classifier here
    
def tissueClass = "Tissue" // usually set to "Tissue"; whatever class ur tissue selector function runs on
def ignoreClass = "Not Lung" // set to class to ignore in tissue selection

double coverage = 0.02 // % total sample area covered by random ROIs(decimal form)

// dims (microns) of ROIs
    double ROI_x_um = 500
    double ROI_y_um = 500

int reps = 3 //how many replicates

//finish param setup; shouldn't have to edit anything from here out

//////////////////////

//next: skip slide if metadata doesn't match
def entryValue = entry.getMetadataValue(specifiedKey)
if ((entryValue != specifiedValue)||(entry.getMetadataValue("Type") != "Sample")) { // also excludes anything that's not a sample; in this case, samples have metadata type=sample
    print "Skipping ${entry.getImageName()}" // comment out if you don't want to see skips
    return 
}

print "Running on ${entry.getImageName()}"

// next: make tissue annotation if one does not already exist
resetSelection()
def containsTissueAnno = getAnnotationObjects().any { it.getPathClass()?.getName() == tissueClass }
if (!containsTissueAnno) {
   print "Creating tissue annotation" 
    createAnnotationsFromPixelClassifier(tissueClassifier, 1000000.0, 5000000.0) // run tissue selector if Tissue anno doesn't already exist
}
print "Tissue annotation already exists for image"

// next: clunky script to remove any areas labeled ignoreClass from tissue selection
selectObjectsByClassification(ignoreClass)
mergeSelectedAnnotations()
def grossTissue = getAnnotationObjects().find{it.getPathClass() == getPathClass(tissueClass)}
def grossTissueGeom = grossTissue.getROI().getGeometry()
def ignore = getAnnotationObjects().find{it.getPathClass() == getPathClass(ignoreClass)}
if (ignore != null){
    def ignoreGeom = ignore.getROI().getGeometry()
    def difference = grossTissueGeom.difference(ignoreGeom)
    selectObjectsByClassification(tissueClass) 
    removeSelectedObjects()
    def tRoi = GeometryTools.geometryToROI(difference, grossTissue.getROI().getImagePlane())
    def tissueanno = PathObjects.createAnnotationObject(tRoi, getPathClass(tissueClass))
    addObject(tissueanno)
    selectObjects(tissueanno)
    }
//

// next: random rectangles script from Krstin Gallik
//setup some classes that will be used later and selecting the main annotation that the random rectangles will be in
def annoList = getAnnotationObjects().findAll{it.getPathClass() == getPathClass(tissueClass)}
def anno = annoList.get(0)
def tissue = anno.getROI()
def tissuePathClass = getPathClass(tissueClass)
def roiClass = getPathClass("ROI")
// I needed to analyze 20% of the area from randomly selected ROIs, this part is to calculate how many 2048x2048 rectangles will cover 20% of the annotation
def tissue_area = tissue.getArea() //this number is in pixels^2
double coverage_area = tissue_area*coverage
double micronsPerPixel = getCurrentImageData().getServer().getPixelCalibration().getAveragedPixelSize()
double ROI_x_px = ROI_x_um/micronsPerPixel // convert width to px
double ROI_y_px = ROI_y_um/micronsPerPixel // ^^ height
double ROI_area_px = ROI_x_px*ROI_y_px
double num_ROIs = Math.round(coverage_area/ROI_area_px)


// start for loop
for (rep in 1..reps) {
anno.each{it.setPathClass(tissuePathClass)} //assigned a classification to the main annotation so it is not included in future steps that are not shown in this code snippet
print "${num_ROIs} ROIs created for replicate ${rep}"
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

// collect all ROIs, run cell detector and classifier
selectObjectsByClassification("ROI") 
runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImage":"DAPI","requestedPixelSizeMicrons":0.5,"backgroundRadiusMicrons":3.0,"backgroundByReconstruction":true,"medianRadiusMicrons":0.0,"sigmaMicrons":2.0,"minAreaMicrons":10.0,"maxAreaMicrons":100.0,"threshold":200.0,"watershedPostProcess":true,"cellExpansionMicrons":3.0,"includeNuclei":true,"smoothBoundaries":true,"makeMeasurements":true}')
runObjectClassifier(cellClassifier)

// merge and end for loop
selectObjectsByClassification("ROI")
mergeSelectedAnnotations()
classifySelected("ROIs")
} // end loop

// names
int counter = 1
for (annotation in getAnnotationObjects().findAll{it.getPathClass() == getPathClass("ROIs")}) {
    annotation.setName("ROI " + counter++)}
    
print("Operation complete for image " + getProjectEntry().getImageName() + "! " + reps + " replicates performed.")

///////////////////////

// troubleshooting
// WARN: No valid (i.e. area or point) annotation to merge
    // probably means you didn't set any annotations of ignoreClass

//////////////////////