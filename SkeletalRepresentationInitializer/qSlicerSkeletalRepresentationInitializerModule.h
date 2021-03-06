/*==============================================================================

  Program: 3D Slicer

  Portions (c) Copyright Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

==============================================================================*/

#ifndef _qSlicerSkeletalRepresentationInitializerModule_h
#define _qSlicerSkeletalRepresentationInitializerModule_h

// Slicer includes
#include "qSlicerLoadableModule.h"

// SkeletalRepresentationInitializer includes
#include "qSlicerSkeletalRepresentationInitializerModuleExport.h"

class qSlicerSkeletalRepresentationInitializerModulePrivate;

/// \ingroup Slicer_QtModules_ExtensionTemplate
class Q_SLICER_QTMODULES_SKELETALREPRESENTATIONINITIALIZER_EXPORT
qSlicerSkeletalRepresentationInitializerModule
  : public qSlicerLoadableModule
{
  Q_OBJECT
#ifdef Slicer_HAVE_QT5
  Q_PLUGIN_METADATA(IID "org.slicer.modules.loadable.qSlicerLoadableModule/1.0");
#endif
  Q_INTERFACES(qSlicerLoadableModule)

public:

  typedef qSlicerLoadableModule Superclass;
  explicit qSlicerSkeletalRepresentationInitializerModule(QObject *parent=nullptr);
  virtual ~qSlicerSkeletalRepresentationInitializerModule() override;

  qSlicerGetTitleMacro(QTMODULE_TITLE)

  virtual QString helpText()const override;
  virtual QString acknowledgementText()const override;
  virtual QStringList contributors()const override;

  virtual QIcon icon()const override;

  virtual QStringList categories()const override;
  virtual QStringList dependencies() const override;

protected:

  /// Initialize the module. Register the volumes reader/writer
  virtual void setup() override;

  /// Create and return the widget representation associated to this module
  virtual qSlicerAbstractModuleRepresentation * createWidgetRepresentation() override;

  /// Create and return the logic associated to this module
  virtual vtkMRMLAbstractLogic* createLogic() override;

protected:
  QScopedPointer<qSlicerSkeletalRepresentationInitializerModulePrivate> d_ptr;

private:
  Q_DECLARE_PRIVATE(qSlicerSkeletalRepresentationInitializerModule)
  Q_DISABLE_COPY(qSlicerSkeletalRepresentationInitializerModule)

};

#endif
