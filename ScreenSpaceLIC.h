#ifndef SCREENSPACEUFLIC_H
#define SCREENSPACEUFLIC_H
#include "UFLIC.h"

template <typename EvalType, typename VecType>
class ScreenSpaceLIC : public UFLIC<EvalType, VecType, 2>
{
public:

  using FieldType = typename UFLIC<EvalType, VecType, 2>::FieldType;

  ScreenSpaceLIC(const vtkm::Id2 &_dim
                 ,int _ttl = 1,
                 int _loop_cnt = 1) :
    UFLIC<EvalType, VecType, 2>(_ttl),
    dim(_dim),
    loop_cnt(_loop_cnt)
  {
    this->propFieldArray[0].Allocate(dim[0] * dim[1]);
    this->propFieldArray[1].Allocate(dim[0] * dim[1]);

    this->omegaArray.Allocate(dim[0] * dim[1]);
    this->texArray.Allocate(dim[0] * dim[1]);

    this->canvasArray.resize(this->ttl);
    for (int i = 0; i < this->ttl; i++) {
      this->canvasArray[i].Allocate(dim[0] * dim[1]);
    }

  }

  template< typename VecFld>
  vtkm::cont::ArrayHandle<FieldType> draw(
                                        std::vector<VecFld> &sl,
                                        std::vector<VecFld> &sr
                                        )
  {
      vtkm::cont::ArrayHandleConstant<vtkm::Id> zero(0, dim[0]*dim[1]);



      for (int i=0; i<this->texArray.GetNumberOfValues(); i++){
        this->texArray.GetPortalControl().Set(i,rand()%255);
      }

      DoNormalize<FieldType> donorm(dim);
      DoSharpen<FieldType> dosharp(dim);
      DoJitter<FieldType> dojitter(dim);
      DrawLineWorklet<FieldType, VecType, 2>  drawline(dim);

      for (int loop = 0; loop < loop_cnt; loop++){
        for (int i=0; i<this->canvasArray[loop%this->ttl].GetNumberOfValues(); i++){
          this->canvasArray[loop%this->ttl].GetPortalControl().Set(i,rand()%255);
        }
        vtkm::cont::ArrayCopy(zero, this->propFieldArray[0]);
        vtkm::cont::ArrayCopy(zero, this->propFieldArray[1]);
        vtkm::cont::ArrayCopy(zero, this->omegaArray);
        for (int i = 0; i < vtkm::Min(this->ttl, loop+1); i++) {
          drawline.Run(this->canvasArray[i],
                       this->propFieldArray[0],
              this->omegaArray, sl[i], sr[i]);
        }

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

      return this->propFieldArray[0];
  }

  vtkm::Id2 dim;
  int loop_cnt;
};
#endif
