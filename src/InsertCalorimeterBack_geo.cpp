//==========================================================================
//  Implementation of backward ECal insert calorimeter
//--------------------------------------------------------------------------
//  Author: Ryan Milton (UCR)
//  Modified by Sebouh Paul (UCR) to make the hole in the same place for each layer
//==========================================================================

#include "DD4hep/DetFactoryHelper.h"
#include <XML/Helper.h>
#include <XML/Layering.h>

using namespace dd4hep;

static Ref_t createDetector(Detector& desc, xml_h handle, SensitiveDetector sens)
{
  xml_det_t     detElem = handle;
  std::string   detName = detElem.nameStr();
  int           detID   = detElem.id();

  xml_dim_t  dim        = detElem.dimensions();
  double     width      = dim.x(); // Size along x-axis
  double     height     = dim.y(); // Size along y-axis 
  double     length     = dim.z(); // Size along z-axis

  xml_dim_t  pos        = detElem.position(); // Position in global coordinates

  Material   air        = desc.material("Air");
  
  // Getting beampipe hole dimensions
  const xml::Component &beampipe_hole_xml = detElem.child(
    _Unicode(beampipe_hole)
  );
  const double hole_radius_1 = dd4hep::getAttrOrDefault<double>(
    beampipe_hole_xml,
    _Unicode(hole_radius_1),
    5.0*cm
  );
  const double hole_radius_2 = dd4hep::getAttrOrDefault<double>(
    beampipe_hole_xml,
    _Unicode(hole_radius_2),
    2.3*cm
  );

  const double hole_x_2 = dd4hep::getAttrOrDefault<double>(
    beampipe_hole_xml,
    _Unicode(hole_x_2),
    8.3*cm
  );

    
  // Defining envelope
  Box envelope(width / 2.0, height / 2.0, length / 2.0);

  /*
    Hole initial x-position in local coordinate system
    x = 0 in local coordinate system is x = HcalEndcapPInsert_xposition in global coordinate system
    Hole starts at x = 0 cm in global system
  */
  //const double hole_x_1 = 0 - pos.x();

  /* 
    Cutting out the hole from the envelope
    Cut out hole layer-by-layer to get exact desired shape
    Subtracting a cone is difficult since the cone would have to be angled towards the negative x direction
  */
  Tube electron_hole(0., hole_radius_1, length / 2.);
  SubtractionSolid envelope_with_hole1(
	envelope,
        electron_hole,
        Position(
		 0.-pos.x(),
          0.,
          0
        )
      );
  envelope = envelope_with_hole1;
  Tube hadron_hole(0., hole_radius_2, length / 2.);
  SubtractionSolid envelope_with_hole2(
        envelope_with_hole1,
        hadron_hole,
        Position(
		 hole_x_2-pos.x(),
          0.,
          0
        )
      );
  envelope = envelope_with_hole2;
  Box connecting_hole(hole_x_2/2,hole_radius_2, length/2.);
  SubtractionSolid envelope_with_hole3(
        envelope_with_hole2,
        connecting_hole,
        Position(
		 hole_x_2/2-pos.x(),
          0.,
          0
        )
      );
  envelope=envelope_with_hole3;
  

  // Defining envelope volume
  Volume envelopeVol(detName, envelope, air); 
  // Setting envelope attributes
  envelopeVol.setAttributes(
    desc,
    detElem.regionStr(),
    detElem.limitsStr(),
    detElem.visStr()
  );

  PlacedVolume pv;

  
  double z_distance_traversed = 0.;
  int layer_num = 1;
  
  // Looping through all the different layer sections
  for(xml_coll_t c(detElem,_U(layer)); c; ++c)
  {
    xml_comp_t x_layer = c;
    int repeat = x_layer.repeat();  
    double layer_thickness = x_layer.thickness();
    
    // Looping through the number of repeated layers in each section
    for(int i = 0; i < repeat; i++)
    {
      std::string layer_name = detName + _toString(layer_num, "_layer%d");

      Box layer(width / 2., height / 2., layer_thickness / 2.);




      // Removing beampipe shape from each layer
      //Tube layer_hole(0., hole_radius, layer_thickness / 2.);
      IntersectionSolid layer_with_hole(
        layer,
        envelope
      );
      Volume layer_vol(layer_name, layer_with_hole, air);

      

      int slice_num = 1;
      double slice_z = -layer_thickness / 2.; // Keeps track of slices' z locations in each layer

      // Looping over each layer's slices
      for(xml_coll_t l(x_layer,_U(slice)); l; ++l) 
      {
        xml_comp_t x_slice = l;
        double slice_thickness = x_slice.thickness();
        std::string slice_name = layer_name + _toString(slice_num, "slice%d");
        Material slice_mat = desc.material(x_slice.materialStr());		
        slice_z += slice_thickness/2.; // Going to slice halfway point

        // Each slice within a layer has the same hole radius
        Box slice(width/2., height/2., slice_thickness/2.);
        
        IntersectionSolid slice_with_hole(
          slice,
          envelope
        );
        Volume slice_vol (slice_name, slice_with_hole, slice_mat);
        
        // Setting appropriate slices as sensitive
        if(x_slice.isSensitive())
        {
          sens.setType("calorimeter");
          slice_vol.setSensitiveDetector(sens);
        }
        
        // Setting slice attributes
        slice_vol.setAttributes(
          desc,
          x_slice.regionStr(),
          x_slice.limitsStr(),
          x_slice.visStr()
        );

        // Placing slice within layer
        pv = layer_vol.placeVolume(
          slice_vol,
          Transform3D(
            RotationZYX(0, 0, 0),
            Position(
              0.,
              0.,
              slice_z
            )
          )
        );
        pv.addPhysVolID("slice", slice_num);
        slice_z += slice_thickness/2.;
        z_distance_traversed += slice_thickness;
        ++slice_num;
      }

      // Setting layer attributes
      layer_vol.setAttributes(
        desc,
        x_layer.regionStr(),
        x_layer.limitsStr(), 
        x_layer.visStr()
      );
      // Placing each layer inside the envelope volume
      // -length/2. is front of detector in global coordinate system
      // + (z_distance_traversed - layer_thickness) goes to the front of each layer
      // + layer_thickness/2. places layer in correct spot
      // Example: After placement of slices in first layer, z_distance_traversed = layer_thickness
      //          Subtracting layer_thickness goes back to the front of the first slice (Now, z = -length/2)
      //          Adding layer_thickness/2. goes to half the first layer thickness (proper place to put layer)
      //          Each loop over repeat will increases z_distance_traversed by layer_thickness
      pv = envelopeVol.placeVolume(
        layer_vol,
        Transform3D(
          RotationZYX(0, 0, 0),
          Position(
            0.,
            0.,
            -length/2.  + (z_distance_traversed - layer_thickness) + layer_thickness/2.
          )
        )
      );
      pv.addPhysVolID("layer", layer_num);
      layer_num++;
      //hole_x_tracker -= 0.0; // The hole shifts along -x by 0 cm every layer
    }
  }

  DetElement   det(detName, detID);  
  Volume motherVol = desc.pickMotherVolume(det);
  
  // Placing insert in world volume
  auto tr = Transform3D(Position(pos.x(), pos.y(), pos.z() + length / 2.));
  PlacedVolume phv = motherVol.placeVolume(envelopeVol, tr);
  phv.addPhysVolID("system", detID);
  det.setPlacement(phv);

  return det;
}
DECLARE_DETELEMENT(InsertCalorimeterBack, createDetector)
