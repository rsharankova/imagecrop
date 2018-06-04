1;95;0c#include "crop_functions.h"
#include "TRandom3.h"

// larlite
#include "LArUtil/Geometry.h"
#include "LArUtil/LArProperties.h"

// function to generate cropped regions
std::vector< larcv::Particle > generate_regions( const int rows, const int cols,
					    const larcv::ImageMeta& sourceimg_meta, const std::vector<larcv::Image2D>& src_v,
					    const int num_regions, const std::vector<float>& min_occupancy_fraction, const std::vector<float>& thresholds,
					    const int maxattempts, const int randseed, std::vector<larcv::Image2D>& previmg_v ) {

  // we generate random positions in the detector
  // we accept only if there is a minimum occupany in the image
  // or a max attempt limit is reached

  std::vector< larcv::Particle > roi_v;
  
  std::vector< std::vector<float> > vertices;

  int seed = randseed;
  if ( randseed<0 )
    seed = int(std::time(NULL));

  // determine limits based on row,col boundaries
  float minz = 0.3*(0.5*cols); 
  float maxz = 1036.0 - minz;
  const larutil::Geometry* geo = larutil::Geometry::GetME();
  const larutil::LArProperties* larp = larutil::LArProperties::GetME();

  //std::cout << "defining crop regions for source image with: (" << sourceimg_meta.rows() << "," << sourceimg_meta.cols() << ")" << std::endl;
  
  TRandom3 rand( seed );
  int numattempts = 0;
  while ( numattempts<maxattempts && roi_v.size()<num_regions ) {

    // generate random position
    Double_t pos[3];
    pos[0] = rand.Uniform()*256.0;
    pos[1] = -117.0 + 2*117.0*rand.Uniform();
    pos[2] = minz + (maxz-minz)*rand.Uniform();

    // define regions in the planes
    std::vector< std::vector<int> > colranges;
    for (int p=0; p<3; p++) {
      float centerwire = geo->WireCoordinate( pos, p );
      float minwire = centerwire-0.5*cols*sourceimg_meta.pixel_width();
      float maxwire = centerwire+0.5*cols*sourceimg_meta.pixel_width();

      if ( maxwire>sourceimg_meta.max_x() ) {
	float dwire = maxwire-sourceimg_meta.max_x()+1;
	maxwire -= dwire;
	minwire -= dwire;
      }
      else if ( minwire<sourceimg_meta.min_x() ) {
	float dwire = sourceimg_meta.min_x()-minwire+1;
	maxwire += dwire;
	minwire += dwire;
      }
      //std::cout << "centerwire " << centerwire <<" minwire "<< minwire <<" maxwire "<< maxwire << std::endl;
      
      std::vector<int> range(2);
      range[0] = sourceimg_meta.col(minwire);
      range[1] = sourceimg_meta.col(maxwire);
      if ( range[1]-range[0]<cols ) range[1] = range[0]+cols;
      colranges.push_back( range );
      //std::cout << "col min " << colranges[p][0] <<" col max "<< colranges[p][1] <<" "<< colranges.size() << std::endl;
    }
    
    // set the time bounds
    std::vector<int> rowrange(2);
    float centertick = 3200.0 + pos[0]/(larp->DriftVelocity()*0.5);
    float mintick = centertick - 0.5*rows*sourceimg_meta.pixel_height();
    float maxtick = centertick + 0.5*rows*sourceimg_meta.pixel_height();
    if ( maxtick > sourceimg_meta.max_y() ) {
      float dtick = maxtick - sourceimg_meta.max_y() + 1;
      maxtick -= dtick;
      mintick -= dtick;
    }
    else if ( mintick < sourceimg_meta.min_y() ) {
      float dtick = sourceimg_meta.min_y() - mintick + 1;
      maxtick += dtick;
      mintick += dtick;
    }
    //std::cout << "centertick " << centertick <<" mintick "<< mintick <<" maxtick "<< maxtick << std::endl;
    
    rowrange[0] = sourceimg_meta.row( mintick );
    rowrange[1] = sourceimg_meta.row( maxtick );
    if ( rowrange[1]-rowrange[0]!=rows )
      rowrange[1] = rowrange[0]+rows;
    //std::cout << "row min " << rowrange[0] <<" "<< sourceimg_meta.row( mintick )<<" row max "<< rowrange[1] << std::endl;

    // generate ranges, now we determine if requirements are satisfied
    int planes_passing = 0;
    int planes_overlap = 0;
    std::vector<float> occfrac(3,0);
    for (int p=0; p<3; p++) {
      const larcv::Image2D& src = src_v[p];
      larcv::Image2D& previmg = previmg_v[p];
      int abovethresh = 0;
      int overlap = 0;
      for (int r=rowrange[0]; r<rowrange[1]; r++) {
	for (int c=colranges[p][0]; c<colranges[p][1]; c++) {
	  if ( src.pixel(r,c)>thresholds[p] )
	    abovethresh++;
	  if (previmg.pixel(r,c)>0.0)
	    overlap++;
	}
      }
      occfrac[p] = float(abovethresh)/float(rows*cols);
      if (  occfrac[p] > min_occupancy_fraction[p] )
	planes_passing++;
      //std::cout <<"plane "<< p << " occupancy frac " << occfrac[p] << planes_passing << std::endl;
      if( float(overlap)/float(rows*cols)>0.5)
	planes_overlap++;
      std::cout <<"overlapping planes: "<< overlap <<" "<< planes_overlap << std::endl;
    } 
    //std::cout <<"overlapping planes: "<< planes_overlap << std::endl;
    std::cout <<"overlapping planes: "<< previmg_v[0].pixel(0,0) << std::endl;

    
    if ( planes_passing==src_v.size() && planes_overlap==0 ) {
      // create bounding boxes and ROI
      larcv::Particle roi;
      for ( int p=0; p<(int)src_v.size(); p++ ) {
	
	larcv::ImageMeta meta( sourceimg_meta.pos_x( colranges[p][0] ), sourceimg_meta.pos_y( rowrange[0] ),
			       sourceimg_meta.pos_x( colranges[p][1] ), sourceimg_meta.pos_y( rowrange[1] ),
			       rows, cols,
			       src_v[p].meta().id(),
			       (larcv::DistanceUnit_t)(larcv::kUnitWireTime));
	//create bb
	roi.boundingbox_2d( meta ,meta.id() );
	//std::cout << meta.dump() << std::endl;
	previmg_v[p].clear();
      }
      
      roi_v.emplace_back( std::move(roi) );
    }
    numattempts++;
    // std::cout << "attempt " << numattempts << " planes_passing=" << planes_passing << std::endl;
    // for (int p=0; p<3; p++)
    //   std::cout << "  occupancy plane " << p << ": " << occfrac[p] << std::endl;
    
  }//end of attempt loop

  return roi_v;
}


void make_cropped_label_image( const std::vector<larcv::Image2D>& srcimgs,
			       const std::vector<larcv::Image2D>& croppedimgs,
			       const std::vector<larcv::Image2D>& piximgs,
			       const std::vector<larcv::Image2D>& visimgs,
			       const std::vector<float> adcthresholds,
			       std::vector<larcv::Image2D>& labelimg_v, std::vector<larcv::Image2D>& matchimg_v,
			       std::vector<larcv::Image2D>& weightimg_v , bool& skip) {
  

  //std::cout << labelimg_v.size() << std::endl;
  //std::cout << labelimg_v[0].meta().dump() << std::endl;
  //std::cout << labelimg_v[0].pixel(0,0) << std::endl;
  
  int targetplanes[3][2] = {{1,2},{0,2},{0,1}}; //for target images order
  // loop over the planes
  for (int p=0; p<3; p++) { 
    skip=false;
    //count visible pixels
    int nonzero[2] = {0,0};
    // get the adc image for the plane. the container stores all three planes. we grab plane p.
    const larcv::Image2D& adcimg    = croppedimgs.at(p);
    const larcv::ImageMeta& adcmeta = adcimg.meta();
    const larcv::ImageMeta& srcmeta = srcimgs.at(p).meta();

    const larcv::ImageMeta* target_meta[2];
    target_meta[0] = &croppedimgs.at( targetplanes[p][0] ).meta();
    target_meta[1] = &croppedimgs.at( targetplanes[p][1] ).meta();

    int count=0;
      
    // loop over rows and columns of the cropped adc image, as it is a subset
    for (size_t radc=0; radc<adcmeta.rows(); radc++) {
      
      float tick = adcmeta.pos_y( radc );
      
      int rid = srcmeta.row(tick);
      
      for (size_t cadc=0; cadc<adcmeta.cols(); cadc++) {
	float wire = adcmeta.pos_x( cadc );
	
	int cid    = srcmeta.col(wire);	  
	
	// get the adc value
	float adc = adcimg.pixel(radc,cadc);
	
	// if below threshold, skip
	//if ( adc<adcthresholds[p] )  continue;

	for(int path=0; path<2; path++){
            
	  //src images where pixels contain pix label & visibility
	  const larcv::Image2D& piximg    = piximgs.at(2*p+path);
	  const larcv::Image2D& visimg    = visimgs.at(targetplanes[p][path]);
	  
	  // get the label pixel
	  float pix = piximg.pixel(rid,cid);
	  float vis = visimg.pixel(rid,cid);
	  float targetwire = wire + pix;

	  if ( targetwire<target_meta[path]->min_x() || targetwire>=target_meta[path]->max_x() ) {
	    // outside the target cropped image
	    labelimg_v[2*p+path].set_pixel( radc, cadc, -3000.0 );
	    matchimg_v[2*p+path].set_pixel( radc, cadc, 0.0 );
	  }
	  else {
	    // inside. need to adjust the flow.
	    int target_col = target_meta[path]->col( targetwire );
	    float target_adc = croppedimgs.at( targetplanes[p][path] ).pixel( radc, target_col );

	    float target_pix = target_col - (float)(cadc);
	    //std::cout << "inside "<< targetwire <<" "<<target_col <<" "<< target_pix<<std::endl;
	    labelimg_v[2*p+path].set_pixel( radc, cadc, target_pix );
	    matchimg_v[2*p+path].set_pixel( radc, cadc, 1.0 );
	    //std::cout << labelimg_v[2*p+path].pixel(radc, cadc) <<" " << matchimg_v[2*p+path].pixel(radc,cadc) <<" fuck you" << std::endl;
	    if ( target_adc<adcthresholds[ targetplanes[p][path] ] )
	      matchimg_v[2*p+path].set_pixel( radc, cadc, 0.0 );

	  }
	  //zero visibility if flow is weird
	  if(labelimg_v[2*p+path].pixel(radc,cadc)<=-3000.0) matchimg_v[2*p+path].set_pixel( radc, cadc, 0.0 );
	  //count visible pixels
	  if(matchimg_v[2*p+path].pixel(radc,cadc)>0) nonzero[path]++;
	  //debug
	  //count++;
	  //if(count> 10) continue;
	  //std::cout << labelimg_v[2*p+path].pixel(radc, cadc) <<" " << matchimg_v[2*p+path].pixel(radc,cadc) <<" fuck you" << std::endl;

	  
	}//end of path loop

      }//end of col loop      
      
    }// end of row loop
    //skip this image if less than 20 visible pixels in y
    //std::cout<< "plane: "<< p <<" nonzero vis: "<< nonzero[0] <<", "<<nonzero[1] << std::endl;
    if(p==2 && (nonzero[0]<50 || nonzero[1]<50)) skip=true;
  }// end of plane loop


}

