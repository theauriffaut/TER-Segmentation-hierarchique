#ifndef GRAPHVERTEX_H
#define GRAPHVERTEX_H

#include <vector>
#include <map>
#include <algorithm>
#include <iostream>
#include <utility>
#include <iomanip>

/**
 * @brief La classe GraphVertex représente un somme d'un graphe non-orienté.
 */
class GraphVertex
{
public:
    /**
     * @brief Constructeur par défaut
     */
    inline GraphVertex ( ){}

    /**
     * @brief Constructeur prenant en paramètre l'id de la face du primal et l'attribue
     * à cet objet GraphVertex représentant son équivalent dans le graphe dual.
     *
     * @param id: Identifiant OpenMesh de la face dans le primal.
     */
    inline GraphVertex( int id ){
        this->m_id = id;
    }

    /**
     * @brief Liste d'adjacence du GraphVertex. C'est une map composé d'une clé de
     * type int (l'identifiant OpenMesh dans le primal de la face correspondant
     * à un sommet voisin de l'objet courant dans le dual) et d'une valeur de type
     * double, le poids associé à l'arête formé par ces deux sommets.
     */
    std::map<int , double> m_adjacencyList;

    /**
     * @brief Getter de m_id
     */
    inline int getId() const { return this->m_id; }

    /**
     * @brief Fonction d'affichage d'un GraphVertex: Affiche l'id de l'objet courant,
     * suivi de ses voisins et des poids associés aux arêtes les reliant.
     */
    inline void displayVertex() {
        std::cout << m_id << ":->";
        for( std::map<int , double>::iterator it = m_adjacencyList.begin() ; it != m_adjacencyList.end() ; ++it ) {
            std::cout << it->first << "=" << it->second << "->";
        }
        std::cout << std::endl;
    }

private:
    /**
     * @brief id de la face du primal correspondant à l'objet courant
     * dans le dual.
     */
    int m_id;
};

#endif // GRAPHVERTEX_H
