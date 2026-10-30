// sets channel brightness parameters by batch
def entry = getProjectEntry()

// Look for a specific metadata column/key, for example: "Stain"
def stainType = entry.getMetadataValue("BatchScanDate")

// Skip if the metadata doesn't match your target
if (BatchScanDate != "20260622") {
    print "Skipping ${entry.getImageName()}"
    return 
}

setChannelDisplayRange('DAPI', 225, 1200)
setChannelDisplayRange('TRITC', 175, 850)
setChannelDisplayRange('Cy5', 238, 825)