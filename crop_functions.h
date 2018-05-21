#include <iostream>
#include <vector>

// larcv
#include "larcv/core/DataFormat/ImageMeta.h"
#include "larcv/core/DataFormat/Image2D.h"
#include "larcv/core/DataFormat/Particle.h"

// larlite
#include "DataFormat/mctrack.h"
#include "DataFormat/mcshower.h"

// function to generate cropped regions
std::vector< larcv::Particle > generate_regions( const int rows, const int cols,
					    const larcv::ImageMeta& sourceimg_meta, const std::vector<larcv::Image2D>& src_v,
					    const int num_regions, const std::vector<float>& min_occupany_fraction, const std::vector<float>& thresholds,
					    const int maxattempts, const int randseed );

void make_cropped_label_image( const std::vector<larcv::Image2D>& origimgs,
			       const std::vector<larcv::Image2D>& croppedimgs, const std::vector<larcv::Image2D>& piximgs, const std::vector<larcv::Image2D>& visimgs,
			       const float adcthreshold,
			       std::vector<larcv::Image2D>& labelimg_v, std::vector<larcv::Image2D>& matchimg_v, std::vector<larcv::Image2D>& weightimg_v );

