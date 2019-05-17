#include "graph.h"

Graph::Graph()
{
}

int Graph::getVertexById( int id ) const {
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
    /**
     * On vérifie l'existence du sommet au préalable et on ne l'ajoute pas si il
     * existe déjà.
     */
    int result = getVertexById( vertexID );
    if( result == -1 ){
        GraphVertex v( vertexID );
        m_vertices.push_back( v );
    }
}

void Graph::addEdge( int vertexID0 , int vertexID1 , double weight ) {
    /**
     * On utilise la fonction addVertex sur les deux ID passés en paramètres. Ils
     * seront créé si ils n'existent pas.
     */
    addVertex( vertexID0 );
    addVertex( vertexID1 );
    //GraphVertex v0 = ;
    //GraphVertex v1 = ;
    m_vertices[getVertexById( vertexID0 )].m_adjacencyList[vertexID1] = weight;
    m_vertices[getVertexById( vertexID1 )].m_adjacencyList[vertexID0] = weight;
    m_nbEdges++;
}

void Graph::displayGraph( ) {

    for( int i = 0 ; i < m_vertices.size() ; ++i ) {
        m_vertices[i].displayVertex();
    }
}

int Graph::getNbVertices() const
{
    return static_cast<int>(m_vertices.size());
}

int Graph::getNbEdges() const
{
    return m_nbEdges;
}

GraphVertex &Graph::operator[](int i)
{
    return m_vertices[i];
}

void Graph::clear()
{
    m_vertices.clear();
}
