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
    int result = getVertexById( vertexID );
    if( result == -1 ){
        GraphVertex v( vertexID );
        m_vertices.push_back( v );
    }
}

void Graph::addEdge( int vertexID0 , int vertexID1 , double weight ) {
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

std::vector<std::vector<double>> Graph::adjacencyMatrix() {
    std::vector<std::vector<double>> matrix;

    for( int i = 0 ; i < m_vertices.size() ; ++i ) {
        for( std::map<int , double>::iterator it = m_adjacencyList.begin() ; it != m_adjacencyList.end() ; ++it ) {
            matrix[i][it->first] = it->second;
        }
    }
}

std::vector<int> Graph::stoerWagner(){
    std::vector<std::vector<double>> matrix( adjacencyMatrix() );

    int graphSize = matrix.size();
    vector<vector<int>> v( graphSize );
    double bestResult = 1000 * 1000;
    std::vector<int> set;

    for( int i = 0 ; i < graphSize ; ++i ) {
        v[i][1] = i;
    }

    vector<double> w( n );
    vector<bool> exist( n , true ) , in_a( n );

    for( int ph = 0 ; ph < n - 1 ; ph++ ) {
        fill( in_a.begin() , in_a.end() , false );
        fill( w.begin() , w.end() , 0.0 );

        for( int it = 0, prev ; it < n - ph ; it++ ) {
            int sel = -1;
            for( int i = 0 ; i < n ; ++i )
                if( exist[i] && !in_a[i] && ( sel == -1 || w[i] > w[sel] ) )
                    sel = i;

            if( it == n - ph - 1) {
                if ( w[sel] < bestResult ) {
                    bestResult = w[sel];
                    set = v[sel];
                }

                v[prev].insert(v[prev].end(), v[sel].begin(), v[sel].end());
                for( int i = 0 ; i < n ; i++)
                    matrix[prev][i] = matrix[i][prev] += matrix[sel][i];

                exist[sel] = false;
            }

            else {
                in_a[sel] = true;
                for( int i = 0 ; i < n ; ++i )
                    w[i] += matrix[sel][i];
                prev = sel;
            }
        }
    }

    return set;
}
