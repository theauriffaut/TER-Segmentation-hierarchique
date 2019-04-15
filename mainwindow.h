#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QFileDialog>
#include <QMainWindow>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <utility>
#include <map>

namespace Ui {
class MainWindow;
}

using namespace OpenMesh;
using namespace OpenMesh::Attributes;
using namespace std;

struct MyTraits : public OpenMesh::DefaultTraits
{
    // use vertex normals and vertex colors
    VertexAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color );
    // store the previous halfedge
    HalfedgeAttributes( OpenMesh::Attributes::PrevHalfedge );
    // use face normals face colors
    FaceAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color );
    EdgeAttributes( OpenMesh::Attributes::Color );
    // vertex thickness
    VertexTraits{float thickness;};
    // edge thickness
    EdgeTraits{float thickness;};
};
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyMesh;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:

    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

    void displayMesh(MyMesh *_mesh);
    void resetAllColorsAndThickness(MyMesh* _mesh);

    float angleFF(MyMesh *_mesh, int faceID0, int faceID1);
    MyMesh::Point faceGravityCenter( MyMesh *_mesh , int faceID );
    float geodesicDistance( MyMesh *_mesh , int faceID0 , int faceID1 );
    void computeGeodesicDistances( MyMesh *_mesh);
    void computeAngularDistances( MyMesh *_mesh );
    void segmentationSimple(MyMesh* _mesh, int k);

    float avgAngularDistance();
    float avgGeodesicDistances();

private slots:
    void on_pushButton_chargement_clicked();
    void on_pushButton_segmentation_clicked();

private:
    MyMesh mesh;
    std::map<pair<int, int>, float> angularDistances;
    std::map<pair<int, int>, float> geodesicDistances;
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
