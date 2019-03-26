#ifndef SCREENSPACEUFLIC_H
#define SCREENSPACEUFLIC_H
#include "UFLIC.h"

template <typename EvalType, typename VecType>
class ScreenSpaceLIC : public UFLIC<EvalType, VecType, 2>
{
public:

  using FieldType = typename UFLIC<EvalType, VecType, 2>::FieldType;
  typedef EulerIntegrator<EvalType, VecType, 2> IntegratorType;
  typedef ParticleAdvectionWorklet<IntegratorType, VecType, 2> ParticleAdvectionWorkletType;

  ScreenSpaceLIC(const vtkm::Id2 &_dim
                 ,vtkm::Float32 _s = 1.0
                 ,int _ttl = 1
                 ,int _loop_cnt = 1) :
    UFLIC<EvalType, VecType, 2>(_ttl),
    dim(_dim)
  ,stepsize(_s)
  ,loop_cnt(_loop_cnt)
  {
    std::cout << stepsize << std::endl;
    this->do_print = true;
    this->propFieldArray[0].Allocate(dim[0] * dim[1]);
    this->propFieldArray[1].Allocate(dim[0] * dim[1]);

    this->omegaArray.Allocate(dim[0] * dim[1]);
    this->texArray.Allocate(dim[0] * dim[1]);

    this->canvasArray.resize(this->ttl);
    for (int i = 0; i < this->ttl; i++) {
      this->canvasArray[i].Allocate(dim[0] * dim[1]);
    }

  }

  template< typename VecFld, typename VecIdx>
  void draw(
          VecFld &vecArray,
          VecIdx &indexArray,
            vtkm::cont::ArrayHandle<vtkm::Float32> &depth
            )
  {
      vtkm::cont::ArrayHandleConstant<vtkm::Id> zero(0, dim[0]*dim[1]);

      std::vector<vtkm::cont::ArrayHandle<vtkm::Vec<VecType, 2>>> sl(this->ttl), sr(this->ttl);
      vtkm::worklet::DispatcherMapField<ResetParticles<VecType, 2>> resetDispatcher(dim[0]);
      for (int i = 0; i < this->ttl; i++) {
          sl[i].Allocate(indexArray.GetNumberOfValues());
          resetDispatcher.Invoke(indexArray, sl[i]);
          sr[i].Allocate(indexArray.GetNumberOfValues());
          resetDispatcher.Invoke(indexArray, sr[i]);
      }


      for (int i=0; i<this->texArray.GetNumberOfValues(); i++){
        this->texArray.GetPortalControl().Set(i,rand()%255);
      }

      DoNormalize<FieldType> donorm(dim);
      DoSharpen<FieldType> dosharp(dim);
      DoJitter<FieldType> dojitter(dim);
      DrawLineWorklet<FieldType, VecType, 2>  drawline(dim, stepsize);
      auto  spacing = vtkm::Vec<VecType,2>(1,1);
      for (int loop = 0; loop < loop_cnt; loop++){
        EvalType eval(0, Bounds(0, dim[0], 0, dim[1]), spacing);
        IntegratorType integrator(eval, 3.0);

        ParticleAdvectionWorkletType advect(integrator);
        for (int i=0; i<this->canvasArray[loop%this->ttl].GetNumberOfValues(); i++){
          this->canvasArray[loop%this->ttl].GetPortalControl().Set(i,rand()%255);
        }
        vtkm::cont::ArrayCopy(zero, this->propFieldArray[0]);
        vtkm::cont::ArrayCopy(zero, this->propFieldArray[1]);
        vtkm::cont::ArrayCopy(zero, this->omegaArray);
        for (size_t i = 0; i < vtkm::Min(this->ttl, loop+1); i++) {
          advect.Run(sl[i], sr[i], vecArray);
          drawline.Run(this->canvasArray[i],
                       this->propFieldArray[0],
              this->omegaArray, depth,indexArray, sl[i], sr[i]);
        }
        sr.swap(sl);

        donorm.Run(this->propFieldArray[0],
            this->omegaArray,
            this->propFieldArray[1]);
        if (this->do_print){
          std::stringstream fn;
          fn << "uflic-" << loop << ".pnm";
         this->saveAs(fn.str().c_str(), this->propFieldArray[1], dim[0], dim[1]);
        }

        //REUSE omegaArray as a temporary cache to sharpen
        dosharp.Run(this->propFieldArray[1],
            this->omegaArray);
        dojitter.Run(this->omegaArray,
                     this->texArray,
                     this->canvasArray[(loop) % this->ttl]);

      }
      this->result = this->propFieldArray[1];
      std::stringstream fn;
      fn << "uflic-final" << ".pnm";
      this->saveAs(fn.str().c_str(), this->propFieldArray[1], dim[0], dim[1]);

  }

  vtkm::Float32 stepsize;
  vtkm::Id2 dim;
  int loop_cnt;
};
#endif
