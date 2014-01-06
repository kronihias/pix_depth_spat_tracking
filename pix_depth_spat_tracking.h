/*-----------------------------------------------------------------
LOG
    GEM - Graphics Environment for Multimedia

    Get pixel information

    Copyright (c) 1997-1998 Mark Danks. mark@danks.org
    Copyright (c) Günther Geiger. geiger@epy.co.at
    Copyright (c) 2001-2011 IOhannes m zmölnig. forum::für::umläute. IEM. zmoelnig@iem.at
    Copyright (c) 2002 James Tittle & Chris Clepper
    For information on usage and redistribution, and for a DISCLAIMER OF ALL
    WARRANTIES, see the file, "GEM.LICENSE.TERMS" in this distribution.

-----------------------------------------------------------------*/

#ifndef _INCLUDE__GEM_PIXES_pix_depth_spat_tracking_H_
#define _INCLUDE__GEM_PIXES_pix_depth_spat_tracking_H_

#include "Base/GemPixDualObj.h"

/*-----------------------------------------------------------------
-------------------------------------------------------------------
CLASS

    pix_depth_spat_tracking

KEYWORDS
    pix
    
DESCRIPTION


-----------------------------------------------------------------*/


#include <math.h>

#include "Eigen/Dense"
using namespace Eigen;

class GEM_EXTERN pix_depth_spat_tracking : public GemPixObj
{
    CPPEXTERN_HEADER(pix_depth_spat_tracking, GemPixObj);

    public:

	    //////////
    	// Constructor
    	pix_depth_spat_tracking(t_floatarg f);
  
      bool debug_print;
  
    protected:
    	
    	//////////
    	// Destructor
    	virtual ~pix_depth_spat_tracking();

    	//////////
        virtual void 	processImage(imageStruct &image);

        //////////
        void			trigger(void);

        //////////
        void			threshMess(t_float thresh);
  
  
    //////////
    void			resolutionMess(t_float res_x, t_float res_y, t_float res_z);
    
    
    //////////
    void			verticesMess(int argc, t_atom *argv);
  
    //////////
    void			inverseMess(); // print the inverse matrix
    void			visualizeMess(t_float value); // cut spat from output pix
  
    void outputMatrix(Matrix3d M);
  
    void computeInverse();
  
  //////////
	// the buffer
	int           xsize, ysize;      // proposed x/y-sizes
	int           m_xsize,  m_ysize;
	int           m_csize;
	t_atom       *m_buffer;
	int           m_bufsize;
	
	int           oldimagex;
	int           oldimagey;
  
  //////////
	// navigation
	float         m_xstep;
	float         m_ystep;
  
  /////////
	// pointer to the image data
	unsigned char *m_data;
  
    /////////
    int m_res_x, m_res_y, m_res_z;
    
    float m_thresh;
    
    Matrix3d T; // transformation matrix B->B'
    
    Vector3d v1; // origin, 1st base vector
    Vector3d v2;
    Vector3d v3;
    Vector3d v4;
  
    int *m_counter;
  
    bool m_visualize;
  
        //////////
        // The color outlet
        t_outlet    	*m_outlet;
  
        t_outlet    	*m_outlet_temp;

	private:

        //////////
        // Static member callbacks
	static void		triggerMessCallback(void *data);
  static void		threshMessCallback(void *data, t_floatarg thresh);
	static void		resolutionMessCallback(void *data, t_floatarg res_x, t_floatarg res_y, t_floatarg res_z);
	static void		verticesMessCallback(void *data, t_symbol *, int argc, t_atom *argv);
	static void		printMessCallback(void *data);
  static void		inverseMessCallback(void *data);
  static void		visualizeMessCallback(void *data, t_floatarg value);
  
};

#endif	// for header file
