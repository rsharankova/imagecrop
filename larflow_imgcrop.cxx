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
  //std::string MC_FILE    = argv[3];  // mcinfo
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
  //nentries = 5;
  
  for ( int i=0; i<nentries; i++) {
    dataco.goto_entry(i,"larcv");

    // get rse
    int run;
    int subrun;
    int event;
    dataco.get_id( run, subrun, event );
    
    //std::cout << "entry " << i << " (" << run << "," << subrun << "," << event << ")" << std::endl;

    // adc images (event-wise)
    larcv::EventImage2D* ev_img_v = (larcv::EventImage2D*)dataco.get_larcv_data( "image2d", "wire" );

    // larflow images
    larcv::EventImage2D* ev_pix_v  = (larcv::EventImage2D*)dataco.get_larcv_data( "image2d", "larflow" );
    larcv::EventImage2D* ev_vis_v  = (larcv::EventImage2D*)dataco.get_larcv_data( "image2d", "visibility" );

    const std::vector<larcv::Image2D>& img_v = ev_img_v->image2d_array();
    
    if ( img_v.size()==0 )
      continue;
    
    std::vector<float> thresholds(3,10.0);
    int target_occupied_pix = 1000;
    float target_occupancy  = float(target_occupied_pix)/float(512*832);
    std::vector<float> occupancy(3,target_occupancy);
    //std::vector<float> occupancy(3,0.01);
    //occupancy[2] = 0.002;

    //generate images for overlap
    std::vector<larcv::Image2D> previmg_v;
    for(int p=0; p<3; p++){
      larcv::Image2D previmg( img_v.at(p).meta() );
      previmg.paint(0.0);
      previmg_v.emplace_back( std::move(previmg) );
    }

    std::vector<larcv::Particle> roi_v = generate_regions( 512, 832, img_v.front().meta(), img_v, 10, occupancy, thresholds, 200, -1 , previmg_v); 

    //std::cout << "Number of ROIs returned: " << roi_v.size() << std::endl;
    
    // crop and save
    int nroi = 0;

    for ( int iroi=0; iroi<roi_v.size(); iroi++) {
      if ( nroi>10 ) break;

      // set the output event container
      larcv::EventImage2D* ev_out    = (larcv::EventImage2D*)outlarcv.get_data( "image2d", "adc" );
      larcv::EventImage2D* ev_label  = (larcv::EventImage2D*)outlarcv.get_data( "image2d", "label" );
      larcv::EventImage2D* ev_match  = (larcv::EventImage2D*)outlarcv.get_data( "image2d", "match" );      
      larcv::EventImage2D* ev_wireid = (larcv::EventImage2D*)outlarcv.get_data( "image2d", "wireid" );

      larcv::Particle& roi = roi_v.at(iroi);

      //std::cout << roi.dump() << std::endl;
      
      // crop the image
      std::vector<larcv::Image2D> croppedimgs;
      for ( auto const& img : img_v ) {

	larcv::Image2D cropped = img.crop( roi.boundingbox_2d( img.meta().id() ) );	
	croppedimgs.emplace_back( std::move(cropped) );
      }
      std::vector<larcv::Image2D> label_v;
      std::vector<larcv::Image2D> match_v;
      std::vector<larcv::Image2D> wireid_v;      
      for(int p=0; p<3; p++){
	for(int i=0; i<2; i++){
	  // make output label image
	  larcv::Image2D label( croppedimgs.at(p).meta() );
	  label.paint(-3000.0);
	  // make output weight image
	  larcv::Image2D wireid( croppedimgs.at(p).meta() );
	  wireid.paint(-1.0);
	  // make visibility images
	  larcv::Image2D match( croppedimgs.at(p).meta() );
	  match.paint(0.0);
	  //
	  label_v.emplace_back( std::move(label) );
	  match_v.emplace_back( std::move(match) );
	  wireid_v.emplace_back( std::move(wireid) );
	}
	
      }
      bool skip=false;
      make_cropped_label_image( ev_img_v->image2d_array(), croppedimgs, ev_pix_v->image2d_array(), ev_vis_v->image2d_array(),
				thresholds, label_v, match_v, wireid_v, skip );
				
      //std::cout << label_v.size() <<" "<<match_v.size() <<" "<< croppedimgs.size()<< " "<< skip << std::endl;
      //check skip
      bool goodcrop = true;
      for ( int i=0; i<6; i++) {
	if ( label_v[i].meta().cols()!=832 || label_v[i].meta().rows()!=512
	     || match_v[i].meta().cols()!=832 || match_v[i].meta().rows()!=512
	     || wireid_v[i].meta().cols()!=832 || wireid_v[i].meta().rows()!=512)
	  goodcrop = false;
      }
      for ( int i=0; i<3; i++) {
	if ( croppedimgs[i].meta().cols()!=832 || croppedimgs[i].meta().rows()!=512)
	  goodcrop = false;
      }
      
      if(skip==true) continue;
      if(goodcrop==false) continue;

      ev_label->emplace( std::move( label_v ) );
      ev_match->emplace( std::move( match_v ) );
      ev_out->emplace( std::move( croppedimgs ));
      ev_wireid->emplace( std::move( wireid_v ) );

      // set the entry number      
      outlarcv.set_id( run, subrun, event*10 + iroi );

      outlarcv.save_entry();

      nroi++;
    }
    
  }
  
  outlarcv.finalize();

  return 0;
}
