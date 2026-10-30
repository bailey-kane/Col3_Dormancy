///////////////////////

// this script:
// 1. removes all annos labeled Region*
// 2. SKIPS any entries that do not match specified metadata
// by Bailey 8/7/2026

//////////////////////

// set up to run over project each time
def specifiedKey = "BatchScanDate" // what key you want to sort on (batch)
def specifiedValue = "20260424" // what value you want to match to
def classToRemove = "Region*" // usually doesnt change

def entry = getProjectEntry()
def entryValue = entry.getMetadataValue(specifiedKey)

// Skip if the metadata doesn't match target
if ((entryValue != specifiedValue)||(entry.getMetadataValue("Type") != "Sample")) {
    print "Skipping ${entry.getImageName()}"
    return 
}

def imageData = entry.readImageData()
def hierarchy = imageData.getHierarchy()

// removes all annos labeled "Region*"
def targetPathClass = getPathClass(classToRemove)
def toRemove = getAnnotationObjects().findAll{it.getPathClass() == targetPathClass}
hierarchy.removeObjects(toRemove, true)

print "Removed ${toRemove.size()} annotation(s) of class ${classToRemove} from image ${entry.getImageName()}"