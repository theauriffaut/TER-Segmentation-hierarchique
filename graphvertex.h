#ifndef GRAPHVERTEX_H
#define GRAPHVERTEX_H

#include <vector>
#include <map>
#include <algorithm>
#include <iostream>

class GraphVertex
{
public:

    inline GraphVertex ( ){}

    inline GraphVertex( int id ){
        this->m_id = id;
    }

    inline int getId() { return this->m_id; }
    inline std::map<int , double> getAdjacencyList () { return this->m_adjacencyList; }
    inline void displayVertex() {

        std::cout << m_id << ":->";
        for( std::map<int , double>::iterator it = m_adjacencyList.begin() ; it != m_adjacencyList.end() ; it++ ) {
            std::cout << it->first << "=" << it->second;
            if( it != m_adjacencyList.end() ) {
                std::cout << "->";
            }
        }
        std::cout << std::endl;
    }

private:
    int m_id;
    std::map<int , double> m_adjacencyList;
};

#endif // GRAPHVERTEX_H
