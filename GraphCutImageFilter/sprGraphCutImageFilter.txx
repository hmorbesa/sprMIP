#ifndef __sprGraphCutImageFilter_txx
#define __sprGraphCutImageFilter_txx

#include "sprGraphCutImageFilter.h"

#include <time.h>

#include "itkObjectFactory.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkBinaryCrossStructuringElement.h"
#include "itkConstantBoundaryCondition.h"
#include "itkConstShapedNeighborhoodIterator.h"
#include "itkShiftScaleImageFilter.h"
#include "itkMaskedImageToHistogramFilter.h"

#include <boost/config.hpp>
#include <iostream>
#include <string>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include <boost/graph/read_dimacs.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/graph_traits.hpp>


namespace spr
{

template< typename TImage, typename TMask>
GraphCutImageFilter<TImage, TMask>::GraphCutImageFilter()
{
  this->SetNumberOfRequiredInputs(2);
  this->verbose = false;
  this->enhance = false;
}

template< typename TImage, typename TMask>
void GraphCutImageFilter<TImage, TMask>::SetInputImage(const TImage* image)
{
  itk::ImageToImageFilter<TImage,TMask>::SetNthInput(0, const_cast<TImage*>(image));
}

template< typename TImage, typename TMask>
void GraphCutImageFilter<TImage, TMask>::SetInputMask(const TMask* mask)
{
  itk::ImageToImageFilter<TMask,TMask>::SetNthInput(1, const_cast<TMask*>(mask));
}

template< typename TImage, typename TMask>
typename TImage::ConstPointer GraphCutImageFilter<TImage, TMask>::GetInputImage()
{
  return static_cast< const TImage * >
         ( this->itk::ProcessObject::GetInput(0) );
}

template< typename TImage, typename TMask>
typename TMask::ConstPointer GraphCutImageFilter<TImage, TMask>::GetInputMask()
{
  return static_cast< const TMask * >
         ( this->itk::ProcessObject::GetInput(1) );
}

template< typename TImage, typename TMask>
void GraphCutImageFilter<TImage, TMask>::GenerateData()
{
  typedef itk::Statistics::MaskedImageToHistogramFilter< TImage, TMask >    MaskedHistogramType;
  typedef typename MaskedHistogramType::HistogramMeasurementVectorType      HistogramMeasurementVectorType;
  typedef typename MaskedHistogramType::HistogramSizeType                   HistogramSizeType;
  typedef typename MaskedHistogramType::HistogramType                       HistogramType;

  //Initialization and allocations
  typename TImage::ConstPointer input = this->GetInputImage();
  typename TMask::ConstPointer mask = this->GetInputMask();
  typename TImage::Pointer output = this->GetOutput();
  output->SetRegions(input->GetLargestPossibleRegion());
  output->Allocate();

  typedef itk::ConstShapedNeighborhoodIterator<TImage> ShapedIteratorType;
  typename ShapedIteratorType::RadiusType radius;  radius.Fill(1);
  typename ShapedIteratorType::ConstIterator offset;

  unsigned int ind=1;
  float v_max, v_min, v_sink=1.0, v_src=0;
  v_min = static_cast<float>(itk::NumericTraits< typename TImage::PixelType >::max());
  v_max = static_cast<float>(itk::NumericTraits< typename TImage::PixelType >::min());

  float v_i,v_j;
  bool isInBounds;
  std::vector<typename TImage::IndexType> idxList;

  long n, head, tail;                    /*  number of edges and nodes */
  float cap;

  vertex_color_t color_object;
  typename graph_traits<Graph>::vertex_descriptor src;
  typename graph_traits<Graph>::vertex_descriptor sink;

  //End of initializations.

  //First loop, for assigning linear index to mask image, get the number of active spels
  itk::ImageRegionIterator<TImage> outputIterator(output, output->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<TImage> inputIterator(input, input->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<TMask> maskIterator(mask, mask->GetLargestPossibleRegion());

  for ( inputIterator.GoToBegin(), maskIterator.GoToBegin(), outputIterator.GoToBegin();
       !outputIterator.IsAtEnd(); ++inputIterator, ++maskIterator, ++outputIterator )
  {
    if( maskIterator.Get()>0 )
    {
      idxList.push_back( maskIterator.GetIndex() );
      outputIterator.Set(ind);
      ind++;
      v_min = std::min( v_min, static_cast<float>(inputIterator.Get()) );
      v_max = std::max( v_max, static_cast<float>(inputIterator.Get()) );
    }
    else
      outputIterator.Set(0);
  }
  n = ind+1;

  if( this->enhance )
  {
    HistogramSizeType histogramSize( 1 );
    histogramSize[0] = v_max-v_min+1;  // number of bins for the channel
    typename MaskedHistogramType::Pointer histogramGenerator = MaskedHistogramType::New();
    histogramGenerator->SetInput( input );
    histogramGenerator->SetMaskImage( mask );
    histogramGenerator->SetMaskValue( 1 );
    histogramGenerator->SetAutoMinimumMaximum( true );
    histogramGenerator->SetHistogramSize(histogramSize);
    histogramGenerator->Update();
    v_max = 0.2*histogramGenerator->GetOutput()->Quantile(0,0.98);
  }

  std::vector<bool> labels(n);

  //Second loop, for building the graph
  ShapedIteratorType neighIt(radius, output, output->GetLargestPossibleRegion());
  neighIt.SetNeedToUseBoundaryCondition( true );
  if (input->ImageDimension==2)
  {
      {
        typename ShapedIteratorType::OffsetType off;
        off[0] = 0; off[1] = 1;
        neighIt.ActivateOffset(off);
      }
      {
        typename ShapedIteratorType::OffsetType off;
        off[0] = 0; off[1] = -1;
        neighIt.ActivateOffset(off);
      }
      {
        typename ShapedIteratorType::OffsetType off;
        off[0] = -1; off[1] = 0;
        neighIt.ActivateOffset(off);
      }
      {
        typename ShapedIteratorType::OffsetType off;
        off[0] = 1; off[1] = 0;
        neighIt.ActivateOffset(off);
      }
  }
  if ( input->ImageDimension==3 )
  {
    {
      typename ShapedIteratorType::OffsetType off;
      off[0] = 0; off[1] = 0; off[2] = 1;
      neighIt.ActivateOffset(off);
    }
    {
      typename ShapedIteratorType::OffsetType off;
      off[0] = 0; off[1] = 0; off[2] = -1;
      neighIt.ActivateOffset(off);
    }
    {
      typename ShapedIteratorType::OffsetType off;
      off[0] = 0; off[1] = -1; off[2] = 0;
      neighIt.ActivateOffset(off);
    }
    {
      typename ShapedIteratorType::OffsetType off;
      off[0] = 0; off[1] = 1; off[2] = 0;
      neighIt.ActivateOffset(off);
    }
    {
      typename ShapedIteratorType::OffsetType off;
      off[0] = -1; off[1] = 0; off[2] = 0;
      neighIt.ActivateOffset(off);
    }
    {
      typename ShapedIteratorType::OffsetType off;
      off[0] = 1; off[1] = 0; off[2] = 0;
      neighIt.ActivateOffset(off);
    }
  }

  if( this->verbose )
    std::cout << "Number of nodes: " << n << std::endl;

  //Define the graph and its properties:
  Graph g;  
  property_map < Graph, edge_capacity_t >::type
          capacity     = get(edge_capacity, g);
  property_map < Graph, edge_reverse_t >::type
          reverse_edge = get(edge_reverse, g);
  //Create graph's vertices:
  std::vector<typename graph_traits<Graph>::vertex_descriptor> verts;
  graph_traits < Graph >::vertex_iterator u_iter, u_end;
  for (long vi = 0; vi < n; ++vi) //number of spels + source + sink
  {
      verts.push_back(boost::add_vertex(g));
  }
  src = 0; sink = ind;

  //Add links for all spels inside the mask
  for ( typename std::vector<typename TImage::IndexType>::iterator idx = idxList.begin();
        idx != idxList.end(); ++idx )
  {
    tail = neighIt.GetCenterPixel();
    inputIterator.SetIndex( *idx );
    //value normalization [0 1]
    v_i = (float(inputIterator.Get()) - v_min)/(v_max-v_min);
    v_i = std::min(std::max(double(v_i),0.0),1.0);
    //Add links from node_i to node_j if they both are each other's neighbor
    for (offset = neighIt.Begin(); !offset.IsAtEnd(); ++offset )
    {
      //Checking if neighbor is in bounds
      head = neighIt.GetPixel(offset.GetNeighborhoodOffset(),isInBounds);
      if ( isInBounds )
      {
        outputIterator.SetIndex( (*idx)-offset.GetNeighborhoodOffset() );
        head = outputIterator.Get();
        //Checking if neighbor is in the mask
        if ( head > 0 )
        {
          inputIterator.SetIndex( (*idx)-offset.GetNeighborhoodOffset() );
          //value normalization [0 1]
          v_j = (float(inputIterator.Get()) - v_min)/(v_max-v_min);
          v_j = std::min(std::max(double(v_j),0.0),1.0);
          cap  =  pow( v_i-v_j, 2.0 );
          if(this->verbose)
          {
            std::cout << *idx << " ";
            std::cout << (*idx)-offset.GetNeighborhoodOffset() << " ";
            std::cout << cap << std::endl;
          }

          //Add edge from node_i to node_j and viceversa
          edge_descriptor e1, e2;
          bool in1, in2;
          boost::tie(e1, in1) = add_edge(verts[tail], verts[head], g);
          boost::tie(e2, in2) = add_edge(verts[head], verts[tail], g);
          //node_i,node_j affinity
          capacity[e1] = cap;
          capacity[e2] = cap;
          reverse_edge[e1] = e2;
          reverse_edge[e2] = e1;
        }
      }
      //std::cout <<"end\n";
    }
  }

  for(unsigned int iter=0; iter<this->MaximumNumberOfIterations; iter++)
  {
    clock_t init, final; init=clock();        
    if( this->verbose )
        std::cout << iter << "\t" << v_src << "\t" << v_sink << "\t";
    double numberOfChanges = 0;

    //for all spels in the mask
    for ( typename std::vector<typename TImage::IndexType>::iterator idx = idxList.begin();
          idx != idxList.end(); ++idx  )
    {
      neighIt.SetLocation( *idx );
      inputIterator.SetIndex( *idx );
      v_i = (float(inputIterator.Get()) - v_min)/(v_max-v_min);
      v_i = std::min(std::max(double(v_i),0.0),1.0);      
    //Write the edges between spel and source
      tail = src;
      head = neighIt.GetCenterPixel();
      cap  =  pow( v_i-v_src, 2.0 );
      {
        edge_descriptor e1, e2;
        bool in1, in2;
        boost::tie(e1, in1) = boost::add_edge(verts[tail], verts[head], g);
        boost::tie(e2, in2) = boost::add_edge(verts[head], verts[tail], g);
        capacity[e1] = cap;
        capacity[e2] = 0;
        reverse_edge[e1] = e2;
        reverse_edge[e2] = e1;        
      }
    //Write the edges between spel and sink
      tail = neighIt.GetCenterPixel();
      head = sink;
      cap  =  pow( v_i-v_sink, 2.0 );
      {
        edge_descriptor e1, e2;
        bool in1, in2;
        boost::tie(e1, in1) = boost::add_edge(verts[tail], verts[head], g);
        boost::tie(e2, in2) = boost::add_edge(verts[head], verts[tail], g);
        capacity[e1] = cap;
        capacity[e2] = 0;
        reverse_edge[e1] = e2;
        reverse_edge[e2] = e1;        
      }      
    }//end spel's for

    float flow = boykov_kolmogorov_max_flow(g ,src, sink);    
    if(this->verbose)
        std::cout << flow << "\t";
    float n_src=0, n_sink=0; v_src = 0; v_sink = 0; int i=0;
    for (boost::tie(u_iter, u_end) = vertices(g); u_iter != u_end; ++u_iter)
    {      
      if (*u_iter != src & *u_iter != sink)
      {
        inputIterator.SetIndex( idxList[*u_iter] );
        v_i = (float(inputIterator.Get()) - v_min)/(v_max-v_min);
        v_i = std::min(std::max(double(v_i),0.0),1.0);
        if ( boost::get(color_object, g ,*u_iter) == boost::black_color )
        {
          if (!labels[i])
            numberOfChanges++;
          v_sink += v_i; n_sink++; labels[i] = true;
        }
        else
        {
          if (labels[i])
            numberOfChanges++;
          v_src += v_i; n_src++; labels[i] = false;
        }
        i++;
      }
    }
    v_src = v_src/n_src;
    v_sink = v_sink/n_sink;
    final=clock()-init;
    if( this->verbose )
        std::cout << numberOfChanges << "\t" << (double)final / ((double)CLOCKS_PER_SEC) << std::endl;
    if ( numberOfChanges <= this->ValueTolerance  )
    {
        std::cout << "LevelSetGraphCut ended because minimum tolerance was achieved\n";
        break;
    }
  }// end iterations for

  for ( unsigned int i=0; i<idxList.size(); i++ )
  {
    outputIterator.SetIndex( idxList[i] );
    if ( labels[i] )
      outputIterator.Set( 127 );
    else
      outputIterator.Set( 255 );
  }
}

}// end namespace

#endif
