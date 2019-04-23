#ifndef GRAPHVERTEX_H
#define GRAPHVERTEX_H

#include <vector>
#include <map>
#include <algorithm>
#include <iostream>
#include <utility>

class GraphVertex
{
public:

    inline GraphVertex ( ){}

    inline GraphVertex( int id ){
        this->m_id = id;
    }

    std::map<int , double> m_adjacencyList;
    inline int getId() { return this->m_id; }
    inline void displayVertex() {
        std::cout << m_id << ":->";
        for( std::map<int , double>::iterator it = m_adjacencyList.begin() ; it != m_adjacencyList.end() ; ++it ) {
            std::cout << it->first << "=" << it->second << "->";
        }
        std::cout << std::endl;
    }

private:
    int m_id;
};

#endif // GRAPHVERTEX_H
