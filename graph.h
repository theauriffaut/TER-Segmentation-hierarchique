#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include "graphvertex.h"

class Graph
{
public:
    Graph();
    void addVertex( int vertexID );
    void addEdge( int vertexID0 , int vertexID1 , double weight );
    void removeVertex ( int vertexID );
    void mergeVertices ( int vertexID0 , int vertexID1 );

    int getVertexById( int id ) const;
    void displayGraph( );

    int getNbVertices() const;
    int getNbEdges() const;

    std::vector<std::vector<double>> adjacencyMatrix();
    std::pair<std::vector<GraphVertex> , double> minCut();
    std::vector<int> stoerWagner();
    //std::vector<int> stoerWagner();

    GraphVertex &operator[](int i);

    void clear();

private:
    std::vector<GraphVertex> m_vertices;
    int m_nbEdges = 0;
};

#endif // GRAPH_H
