#include <iostream>
#include <string>

// larcv
#include "larcv/core/DataFormat/EventImage2D.h"
#include "larcv/core/DataFormat/Image2D.h"
#include "larcv/core/DataFormat/Particle.h"

// llcv
#include "Base/DataCoordinator.h"

// opencv
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>

// larcv opencv utils
//#include "CVUtil/CVUtil.h"

#include "crop_functions.h"

int main( int nargs, char** argv ) {

  std::cout << "[ Select Tune Sample ]" << std::endl;

  std::string SUP_FILE   = argv[1];  // supera
  std::string OUT_FILE   = argv[2];

  // holds images and mcinfo trees per event
  larlitecv::DataCoordinator dataco;
  dataco.add_inputfile( SUP_FILE, "larcv" );  
  dataco.initialize();

  // output for larcv information
  larcv::IOManager outlarcv( larcv::IOManager::kWRITE );
  outlarcv.set_out_file( OUT_FILE );
  outlarcv.initialize();
  
  int nentries = dataco.get_nentries( "larcv" );
  nentries = 5;
  
  for ( int i=0; i<nentries; i++) {
    dataco.goto_entry(i,"larcv");

    // get rse
    int run;
    int subrun;
    int event;
    dataco.get_id( run, subrun, event );
    
    std::cout << "entry " << i << " (" << run << "," << subrun << "," << event << ")" << std::endl;

    // adc images (event-wise)
    larcv::EventImage2D* ev_img_v = (larcv::EventImage2D*)dataco.get_larcv_data( "image2d", "wire" );

    const std::vector<larcv::Image2D>& img_v = ev_img_v->image2d_array();
    
    if ( img_v.size()==0 )
      continue;
    
    std::vector<float> thresholds(3,10.0);
    std::vector<float> occupancy(3,0.01);
    occupancy[2] = 0.002;
    std::vector<larcv::Image2D> previmg_v;
    for(int p=0; p<3; p++){
      larcv::Image2D previmg( img_v.at(p).meta() );
      previmg.paint(0.0);
      previmg_v.emplace_back( std::move(previmg) );
    }
    
    std::vector<larcv::Particle> roi_v = generate_regions( 512, 512, img_v.front().meta(), img_v, 10, occupancy, thresholds, 100, -1, previmg_v );

    std::cout << "Number of ROIs returned: " << roi_v.size() << std::endl;

    // crop and save
    int nroi = roi_v.size();
    if ( nroi>10 )
      nroi = 10;
    for ( int iroi=0; iroi<nroi; iroi++) {

      // set the output event container
      larcv::EventImage2D* ev_out    = (larcv::EventImage2D*)outlarcv.get_data( "image2d", "adc" );

      larcv::Particle& roi = roi_v.at(iroi);

      //std::cout << roi.dump() << std::endl;
      //std::cout << "pixelwidth: " << roi.BB(0).pixel_width() << " " << roi.BB(0).pixel_height() << std::endl;
      
      // crop the image
      for ( auto const& img : img_v ) {

	larcv::Image2D cropped = img.crop( roi.boundingbox_2d( img.meta().id() ) );	
	ev_out->emplace( std::move(cropped) );
      }
      
      // set the entry number
      outlarcv.set_id( run, subrun, event*10 + iroi );

      outlarcv.save_entry();
      
    }
    
  }
  
  outlarcv.finalize();

  return 0;
}
