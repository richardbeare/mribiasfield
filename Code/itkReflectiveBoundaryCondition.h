/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkReflectiveBoundaryCondition.h,v $
  Language:  C++
  Date:      $Date: 2009-02-19 19:41:23 $
  Version:   $Revision: 1.6 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkReflectiveBoundaryCondition_h
#define __itkReflectiveBoundaryCondition_h
#include "itkNeighborhood.h"
#include "itkImageBoundaryCondition.h"

namespace itk
{

/** \class ReflectiveBoundaryCondition
 * \brief
 * A function object that determines values outside of image boundaries
 * according to periodic (wrap-around) conditions.
 *
 * The input to this function object is a neighborhood iterator.  This boundary
 * condition object is designed to be given as a template argument to a
 * NeighborhoodIterator or any of the NeighborhoodIterator subclasses.
 * 
 * \ingroup DataRepresentation
 * \ingroup ImageObjects
 */
template<class TImage>
class ITK_EXPORT  ReflectiveBoundaryCondition
  : public ImageBoundaryCondition<TImage>
{
public:
  /** Standard class typedefs. */ 
  typedef ReflectiveBoundaryCondition      Self;
  typedef ImageBoundaryCondition<TImage> Superclass;
  
  /** Extract information from the image type. */
  typedef typename Superclass::PixelType        PixelType;
  typedef typename Superclass::PixelPointerType PixelPointerType;
  typedef typename Superclass::IndexType        IndexType;
  typedef typename Superclass::OffsetType       OffsetType;
  typedef typename Superclass::NeighborhoodType NeighborhoodType;
  
  typedef typename Superclass::NeighborhoodAccessorFunctorType 
                                 NeighborhoodAccessorFunctorType;

  /** Extract information from the image type. */
  itkStaticConstMacro(ImageDimension, unsigned int,Superclass::ImageDimension);
  
  /** Default constructor. */
  ReflectiveBoundaryCondition() 
    {
    this->SetReflectionType(3);
    }

  /** Set/Get the type of reflection - 0 =>  c b a a (a-b) (a-c), 1 => c b a
  (a-b) (a-c) , 2=> c b a a b c, 3=> c b a b c
  */
  void SetReflectionType(unsigned K)
    {
    m_ReflectionType = K;
    
    switch (m_ReflectionType)
      {
      case 0:
	m_ReflectionOffset = 0;
	m_HorizontalReflection=true;
	break;
      case 1:
	m_ReflectionOffset = 1;
	m_HorizontalReflection=true;
	break;
      case 2:
	m_ReflectionOffset = 0;
	m_HorizontalReflection=false;
	break;
      case 3:
	m_ReflectionOffset = 1;
	m_HorizontalReflection=false;
	break;
      default:
	// same as 3
	m_ReflectionOffset = 1;
	m_HorizontalReflection=false;
      }
    }

  /** Computes and returns a neighborhood of appropriate values from
   * neighborhood iterator data.. */
  virtual PixelType operator()(const OffsetType& point_index,
                               const OffsetType& boundary_offset,
                               const NeighborhoodType *data) const; 

  /** Computes and returns the appropriate pixel value from
   * neighborhood iterator data, using the functor. */
  virtual PixelType operator()(
      const OffsetType& point_index,
      const OffsetType& boundary_offset,
      const NeighborhoodType *data,
      const NeighborhoodAccessorFunctorType &neighborhoodAccessorFunctor) const;
private:
  unsigned m_ReflectionType;
  unsigned m_ReflectionOffset;  // 0 or 1
  bool m_HorizontalReflection;
};

} // end namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_ReflectiveBoundaryCondition(_, EXPORT, x, y) namespace itk { \
  _(1(class EXPORT ReflectiveBoundaryCondition< ITK_TEMPLATE_1 x >)) \
  namespace Templates { typedef ReflectiveBoundaryCondition< ITK_TEMPLATE_1 x > \
                                         ReflectiveBoundaryCondition##y; } \
  }

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkReflectiveBoundaryCondition+-.h"
#endif

#if ITK_TEMPLATE_TXX
# include "itkReflectiveBoundaryCondition.txx"
#endif

#endif
