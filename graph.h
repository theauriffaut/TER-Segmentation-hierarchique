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
    int getVertexById( int id );
    void displayGraph( );

private:
    std::vector<GraphVertex> m_vertices;
};

#endif // GRAPH_H
