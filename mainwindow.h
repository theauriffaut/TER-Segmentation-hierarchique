#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QFileDialog>
#include <QMainWindow>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <utility>
#include <map>
#include <algorithm>
#include <cmath>
#include "graph.h"

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

    double angleFF(MyMesh *_mesh, int faceID0, int faceID1);
    MyMesh::Point faceGravityCenter( MyMesh *_mesh , int faceID );
    std::vector<MyMesh::Point> discretizeEdge( MyMesh *_mesh , int edgeID );
    double geodesicDistance( MyMesh *_mesh , int edgeID);
    void computeGeodesicDistances( MyMesh *_mesh);
    void computeAngularDistances( MyMesh *_mesh );
    void computeWeight( MyMesh *_mesh , double coefGeod, int patch );
    void segmentationSimple(MyMesh* _mesh, int k);

    double avgAngularDistances();
    double avgGeodesicDistances();

    void displayGeodesicDistances();
    void displayAngularDistances();

    double dijkstraDual(int v1, int v2);
    void dijkstraByREP(MyMesh *_mesh, int IdFaceREP, int k);

    FPropHandleT<QVector<double>> dist;
    FPropHandleT<int> patchId;



    void computeProb();
    void computeDirectDistances(MyMesh *_mesh, int patch);
    void displayDirectDistances();
private slots:
    void on_pushButton_chargement_clicked();
    void on_pushButton_segmentation_clicked();

private:
    MyMesh mesh;
    QVector<QVector<int>> patches;
    QVector<int> ambiguousFaces;
    unsigned int currentId = 0;
    std::map<pair<int, int>, double> angularDistances;
    std::map<pair<int, int>, double> geodesicDistances;
    std::map<pair<int, int>, double> directDistances;
    QVector<MyMesh::Color> colors;
    Graph dual;
    Ui::MainWindow *ui;
    int nbDijkstraDone = 0;
};

#endif // MAINWINDOW_H
