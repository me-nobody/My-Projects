import qupath.lib.objects.PathCellObject
import qupath.lib.objects.PathDetectionObject;
import qupath.lib.objects.hierarchy.PathObjectHierarchy;
import qupath.lib.roi.ROIs;
import qupath.lib.objects.PathAnnotationObject;
import java.util.Collection;
import qupath.lib.gui.tools.MeasurementExporter
import qupath.lib.objects.PathAnnotationObject
import qupath.lib.objects.TMACoreObject
import qupath.lib.objects.PathRootObject
import java.util.Collections
import java.util.HashMap


// Directory to save the CSVs
def outputDir = buildFilePath(PROJECT_BASE_DIR, "Core_Measurements")
mkdirs(outputDir)


// Get the current image data
def imageData = getCurrentImageData()
// Get the hierarchy of objects within the image
def hierarchy = imageData.getHierarchy()


// Get a list of all detected TMA core objects
def tmaCores = getTMACoreList()
int counter = 1
// Iterate through each TMA core in the list

for (core in tmaCores) {    
    coreName = core.getName()
    // Set up export path
    def filePath = buildFilePath(outputDir, "${coreName}.csv")
    def file = new File(filePath)
    headers = "id,name,class,Nucleus.area,Nucleus.perimeter,Nucleus.circularity,Nucleus.Max_caliper,Nucleus.Min_caliper,Nucleus.eccentricity,Nucleus.Hematoxylin.OD.mean,Nucleus.Hematoxyliin.OD.sum,\
           Nucleus.Hematoxylin.OD.std,Nucleus.Hematoxylin.OD.max,Nucleus.Hematoxylin.OD.min,Nucleus.Hematoxylin.OD.range,Nucleus.DAB.OD.mean,Nucleus.DAB.OD.sum,Nucleus.DAB.OD.std,Nucleus.DAB.OD.max,\
           Nucleus.DAB.OD.min,Nucleus.DAB.OD.range,Cell.area,Cell.perimeter,Cell.circularity,Cell.Max_caliper,Cell.Min_caliper,Cell.eccentricity,Cell.Hematoxylin.OD.mean,Cell.Hematoxylin.OD.std,Cell.Hematoxylin.OD.max,\
           Cell.Hematoxylin.OD.min,Cell.DAB.OD.mean,Cell.DAB.OD.std,Cell.DAB.OD.max,Cell.DAB.OD.min,Cytoplasm.Hematoxylin.OD.mean,Cytoplasm.Hematoxylin.OD.std,Cytoplasm.Hematoxylin.OD.max,Cytoplasm.Hematoxylin.OD.min,\
           Cytoplasm.DAB.OD.mean,Cytoplasm.DAB.OD.std,Cytoplasm.DAB.OD.max,Cytoplasm.DAB.OD.min,Nucleus/Cell.area.ratio"
     
     file.append(headers)
     file.append("\n")
    
    // Clear any previous selection
    hierarchy.getSelectionModel().clearSelection()
    
    selectObjects { p -> p == core }
    
    def detectionsInCore = getDetectionObjects().findAll {
        core.getROI().contains(it.getROI().getCentroidX(), it.getROI().getCentroidY())
    }
    
    if (detectionsInCore.isEmpty()) {
        println "Skipping ${coreName} (no detections)"
        counter++
        continue
    }
    
     // Print the name of the current core (optional)
     println 'Processing TMA Core: ' + coreName    
       def detection = detectionsInCore[0]

       for(cell in detectionsInCore) {
            def id = cell.getID().toString()
            file.append(id)
            file.append(",")
            cell.setName(coreName)            
            def name = cell.getName().toString()
            file.append(name)
            file.append(",")
            def cell_class = cell.getClassification()  
            file.append(cell_class)  
            file.append(",")
            def cell_measures = cell.getMeasurements()   
            def cell_m_val = cell_measures.values()
             // Iterate and print the values
             cell_m_val.each { value ->
             file.append(",")
             file.append(value)  
             }  
             file.append("\n")

             println "                          "       
             def id_ = detection.getID()
                 println id_
             detection.setName(coreName)
             def name_ = detection.getName()
                 println name_
             def d_class = detection.getClassification()
                 println d_class         
             def d_measures = detection.getMeasurements()
                 println d_measures
             // def d_m_val = d_measures.values()     
             // Iterate and print the values
             //     d_m_val.each { value ->
             //     println value                     
             //       }
             println "                                "
  
 
    println "Exported measurements for ${coreName} (${detectionsInCore.size()} detections)"
    counter++

   
}
}
