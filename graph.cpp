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

    std::vector<std::vector<double>> matrix ( getNbVertices() , std::vector<double> ( getNbVertices() , 0.0 ) );

    for( int i = 0 ; i < getNbVertices() ; ++i ) {
        for( std::map<int , double>::iterator it = m_vertices[i].m_adjacencyList.begin() ; it != m_vertices[i].m_adjacencyList.end() ; ++it ) {
            matrix[i][getVertexById(it->first)] = it->second;
        }
    }

    return matrix;
}


std::pair<std::vector<GraphVertex> , double > Graph::minCut() {
    std::vector<GraphVertex> A = { m_vertices[0] };

    while( A.size() != getNbVertices() ) {
        double max = 0.0;
        int idMax = 0;

        for( int i = 0 ; i < A.size() ; ++i ) {
            for( std::map<int , double>::iterator it = A[i].m_adjacencyList.begin() ;
                 it != A[i].m_adjacencyList.end() ;
                 ++it ) {

                bool isInA = false;
                for ( int j = 0 ; j < A.size() ; ++j ){
                    if (A[j].getId() == it->first) {
                        isInA = true;
                        break;
                    }
                }

                if ( isInA ) continue;

                if ( it->second > max ) {
                    max = it->second;
                    idMax = getVertexById( it->first );
                }
            }
        }
        A.push_back( m_vertices[idMax] );
    }

    GraphVertex last = A[A.size() - 1];
    GraphVertex nextToLast = A[A.size() - 2];
    double weight = 0.0;

    for( std::map<int , double>::iterator it = last.m_adjacencyList.begin() ;
         it != last.m_adjacencyList.end () ;
         ++it )
        weight += it->second;

    std::pair<std::vector<GraphVertex> , double> ret;
    std::vector<GraphVertex> first = { nextToLast , last };
    ret = std::make_pair( first , weight );
    return ret;
}

void Graph::removeVertex ( int vertexID ) {
    int id = getVertexById(vertexID);
    GraphVertex v = m_vertices[id];

    for( std::map<int , double>::iterator it = v.m_adjacencyList.begin() ;
         it != v.m_adjacencyList.end() ;
         ++it  ) {

        int currentId = getVertexById(it->first);
        std::map<int, double>::iterator toDelete = m_vertices[currentId].getAdjacencyList().find(v.getId());
        m_vertices[currentId].getAdjacencyList().erase(toDelete);
    }

    m_vertices.erase( m_vertices.begin() + id );
}


void Graph::mergeVertices ( int vertexID0 , int vertexID1 ) {
    GraphVertex v0 = m_vertices[getVertexById(vertexID0)];
    GraphVertex v1 = m_vertices[getVertexById(vertexID1)];

    m_vertices[getVertexById(vertexID0)].getMergedVertices().push_back(vertexID1);
    for ( std::vector<int>::iterator it = v1.mergedVertices.begin() ;
          it != v1.mergedVertices.end() ;
          ++it ) {
        m_vertices[getVertexById(vertexID0)].getMergedVertices().push_back(*it);
    }

    for ( std::map<int , double>::iterator it = v1.m_adjacencyList.begin() ;
          it != v1.m_adjacencyList.end() ;
          ++it ) {

        std::map<int , double>::iterator it2 = v0.m_adjacencyList.find(it->first);

        if(it2 != v0.m_adjacencyList.end() ) {
            it2->second += it->second;
        }

        else if( it->first != vertexID0 ){
            addEdge(vertexID0 , it->first , it->second);
        }
    }
    GraphVertex v0COPY = m_vertices[getVertexById(vertexID0)];
    GraphVertex v1COPY = m_vertices[getVertexById(vertexID1)];
    removeVertex( vertexID1 );
}

std::vector<int> Graph::stoerWagner() {
    double bestResult = 1000 * 1000 * 1000;
    std::vector<int> bestCut;
    while (getNbVertices() > 1) {
        std::cout << "Nombre de sommets restants: " << getNbVertices() << std::endl;
        std::pair<std::vector<GraphVertex> , double > currentCut = minCut();
        std::cout << currentCut.first.size();
        std::vector<GraphVertex> vertices = currentCut.first;
        double weight = currentCut.second;

        if ( weight < bestResult ) {
            bestResult = weight;
            bestCut.clear();
            for ( std::vector<GraphVertex>::iterator it = vertices.begin() ;
                  it != vertices.end() ;
                  ++it ) {
                GraphVertex current = *it;
                for ( std::vector<int>::iterator it2 = current.mergedVertices.begin() ;
                      it2 != current.mergedVertices.end() ;
                      ++it2 ) {
                    bestCut.push_back(*it2);
                }
                bestCut.push_back(current.getId());
            }
        }
        mergeVertices( vertices[0].getId() , vertices[1].getId() );
    }
    return bestCut;
}
