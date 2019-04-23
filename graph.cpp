#include "graph.h"

Graph::Graph()
{
}

int Graph::getVertexById( int id ){
    int result = -1;
    for( int i = 0; i < m_vertices.size() ; ++i ){
        if( m_vertices[i].getId() == id ) {
            result = i;
            break;
        }
    }

    return result;
}

void Graph::addVertex( int vertexID ) {
    int result = getVertexById( vertexID );
    if( result == -1 ){
        GraphVertex v( vertexID );
        m_vertices.push_back( v );
    }
}

void Graph::addEdge( int vertexID0 , int vertexID1 , double weight ) {
    addVertex( vertexID0 );
    addVertex( vertexID1 );
    GraphVertex v0 = getVertexById( vertexID0 );
    GraphVertex v1 = getVertexById( vertexID1 );
    v0.getAdjacencyList()[vertexID1] = weight;
    v1.getAdjacencyList()[vertexID0] = weight;
}

void Graph::displayGraph( ) {
    for( int i = 0 ; i < m_vertices.size() ; ++i ) {
        m_vertices[i].displayVertex();
    }
}
