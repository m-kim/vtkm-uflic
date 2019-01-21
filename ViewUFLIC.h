#ifndef VIEWUFLIC_H
#define VIEWUFLIC_H

class VTKM_RENDERING_EXPORT ViewUFLIC : public vtkm::rendering::View3D
{
public:
  ViewUFLIC(const vtkm::rendering::Scene& scene,
         const vtkm::rendering::Mapper& mapper,
         const vtkm::rendering::Canvas& canvas,
         const vtkm::rendering::Camera& camera,
         const vtkm::rendering::Color& backgroundColor = vtkm::rendering::Color(0, 0, 0, 1),
         const vtkm::rendering::Color& foregroundColor = vtkm::rendering::Color(1, 1, 1, 1)) :
         View3D(scene, mapper, canvas, camera, backgroundColor, foregroundColor){}

    void Paint() override
    {
      this->GetCanvas().Activate();
      this->GetCanvas().Clear();

      this->SetupForWorldSpace();
      this->GetScene().Render(this->GetMapper(), this->GetCanvas(), this->GetCamera());

      this->SetupForScreenSpace();

      this->GetCanvas().Finish();
    }

};
#endif
