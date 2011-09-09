/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkReflectiveBoundaryCondition.txx,v $
  Language:  C++
  Date:      $Date: 2007-02-02 14:46:33 $
  Version:   $Revision: 1.8 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkReflectiveBoundaryCondition_txx
#define __itkReflectiveBoundaryCondition_txx
#include "itkConstNeighborhoodIterator.h"
#include "itkReflectiveBoundaryCondition.h"
namespace itk
{
template<class TImage>
typename ReflectiveBoundaryCondition<TImage>::PixelType
ReflectiveBoundaryCondition<TImage>
::operator()(const OffsetType& point_index, const OffsetType& boundary_offset,
             const NeighborhoodType *data) const
{
  typedef typename OffsetType::OffsetValueType OffsetValueType;
  const ConstNeighborhoodIterator<TImage> * iterator
    = dynamic_cast<const ConstNeighborhoodIterator<TImage> *>(data);
  typename TImage::PixelType *ptr, *edgeptr;
  int linear_index = 0;
  unsigned int i;
  
  // Find the pointer of the closest boundary pixel

  // Return the value of the pixel at the closest boundary point.
  for (i = 0; i < ImageDimension; ++i)
    {
    linear_index += (point_index[i] + boundary_offset[i]) * data->GetStride(i);
    }
  ptr = data->operator[](linear_index);
  edgeptr=ptr;
  // Wrap the pointer around the image in the necessary dimensions.  If we have
  // reached this point, we can assume that we are on the edge of the BUFFERED
  // region of the image.  Boundary conditions are only invoked if touching the
  // actual memory boundary.

  // These are the step sizes for increments in each dimension of the image.
  const typename TImage::OffsetValueType * offset_table
    = iterator->GetImagePointer()->GetOffsetTable();
    
  
  for (i = 0; i < ImageDimension; ++i)
    {
    if (boundary_offset[i] != 0)
      { // If the neighborhood overlaps on the low edge, then search
	// in the positive direction
      ptr += (boundary_offset[i] - 1 + m_ReflectionOffset) * offset_table[i];
      }
    }

  if (m_HorizontalReflection)
    {
    return 2*(*edgeptr) - (*ptr);
    }
  else
    {
    return *ptr;
    }
}

template<class TImage>
typename ReflectiveBoundaryCondition<TImage>::PixelType
ReflectiveBoundaryCondition<TImage>
::operator()(const OffsetType& point_index, const OffsetType& boundary_offset,
             const NeighborhoodType *data,
             const NeighborhoodAccessorFunctorType &neighborhoodAccessorFunctor) const
{
  typedef typename OffsetType::OffsetValueType OffsetValueType;
  const ConstNeighborhoodIterator<TImage> * iterator
    = dynamic_cast<const ConstNeighborhoodIterator<TImage> *>(data);
  typename TImage::PixelType *ptr, *edgeptr;
  int linear_index = 0;
  unsigned int i;
  
  // Find the pointer of the closest boundary pixel
  //  std::cout << "Boundary offset = " << boundary_offset << std::endl;
  // std::cout << "point index = " << point_index << std::endl;


  // Return the value of the pixel at the closest boundary point.
  for (i = 0; i < ImageDimension; ++i)
    {
    linear_index += (point_index[i] + boundary_offset[i]) * data->GetStride(i);
    }
  ptr = data->operator[](linear_index);
  edgeptr=ptr;
  // Wrap the pointer around the image in the necessary dimensions.  If we have
  // reached this point, we can assume that we are on the edge of the BUFFERED
  // region of the image.  Boundary conditions are only invoked if touching the
  // actual memory boundary.

  // These are the step sizes for increments in each dimension of the image.
  const typename TImage::OffsetValueType * offset_table
    = iterator->GetImagePointer()->GetOffsetTable();
    
  
  for (i = 0; i < ImageDimension; ++i)
    {
    if (boundary_offset[i] != 0)
      { // If the neighborhood overlaps on the low edge, then wrap from the
      // low edge of the image.
      ptr += (boundary_offset[i] - 1 + m_ReflectionOffset) * offset_table[i];

      }
    }
  if (m_HorizontalReflection)
    {
    return 2*neighborhoodAccessorFunctor.Get(edgeptr) - neighborhoodAccessorFunctor.Get(ptr);
    }
  else 
    {
    return(neighborhoodAccessorFunctor.Get(ptr));
    }
}
} // end namespace itk

#endif
