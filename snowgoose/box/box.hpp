#ifndef __BOX_HPP__
#define __BOX_HPP__

#include <common/smartptr.hpp>
/**
  * Class for handling boxes
  */

namespace snowgoose {
    
template <class FT> struct Box {
  Box(int dim = 0) : mA(dim), mB(dim)
  {
    mDim = dim;
  }

 
  /**
    * Box number of dimensions
    */
  int mDim;

  /**
   * "Left" box bounds
   */
  SmartArrayPtr<FT> mA;

  /**
   * "Right" box bounds
   */
  SmartArrayPtr<FT> mB;
};

}

#endif
