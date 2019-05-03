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

    int getVertexById( int id ) const;
    void displayGraph( );

    int getNbVertices() const;
    int getNbEdges() const;

    GraphVertex &operator[](int i);

    void clear();

private:
    std::vector<GraphVertex> m_vertices;
    int m_nbEdges = 0;
};

#endif // GRAPH_H
