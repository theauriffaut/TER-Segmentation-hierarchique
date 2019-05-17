#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include "graphvertex.h"

/**
 * @brief La classe Graph représente un graphe non-orienté, contenant de multiples
 * objets GraphVertex
 */
class Graph
{
public:
    /**
     * @brief Constructeur par défaut de la classe Graph
     */
    Graph();

    /**
     * @brief Fonction d'ajout d'un vertex au graphe.
     *
     * @param vertexID: Id de la face du primal correspondant au sommet à insérer dans le graphe dual
     */
    void addVertex( int vertexID );

    /**
     * @brief Fonction ajoutant un arête entre deux sommets du graphe. Vérifie au
     * préalable l'existence des sommets et les créé si ils n'existent pas. Ajoute
     * chaque sommet à la liste d'adjacence de l'autre avec le poids associé.
     *
     * @param vertexID0: Id de la face du primal correspondant au premier sommet à relier.
     *
     * @param vertexID1: Id de la face du primal correspondant au second sommet à relier.
     *
     * @param weight: Poids associé à l'arête reliant les deux sommets.
     */
    void addEdge( int vertexID0 , int vertexID1 , double weight );

    /**
     * @brief Renvoie la position dans le tableau m_vertices du sommet dont l'attribut
     * m_id vaut le paramètre id.
     *
     * @param id: Id du sommet à rechercher.
     *
     * @return La position dans le tableau m_vertices du sommet recherché. -1 si ce sommet
     * n'existe pas.
     */
    int getVertexById( int id ) const;

    /**
     * @brief Fonction d'affichage d'un graphe. Itère les objets GraphVertex stockés
     * dans m_vertices et appelle la fonction displayVertex() de chacun d'entre eux.
     */
    void displayGraph( );

    /**
     * @brief Fonction renvoyant la taille du graphe (soit la taille de m_vertices)
     *
     * @return La taille du tableau m_vertices
     */
    int getNbVertices() const;

    /**
     * @brief Getter de l'attribut m_nbEdges
     *
     * @return Valeur de m_nbEdges
     */
    int getNbEdges() const;

    /**
     * @brief Surcharge de l'opérateur [] pour accéder aux éléments du tableau
     * m_vertices avec la syntaxe Graph[i];
     *
     * @param i: Position de l'élément auquel accéder dans m_vertices
     *
     * @return Référence vers l'objet GraphVertex stocké à la position i du tableau m_vertices
     */
    GraphVertex &operator[](int i);

    /**
     * @brief Vide le graphe de ses sommets
     */
    void clear();

private:
    /**
     * @brief Tableau contenant les GraphVertex composant le graphe
     */
    std::vector<GraphVertex> m_vertices;

    /**
     * @brief Nombre d'arêtes du graphe
     */
    int m_nbEdges = 0;
};

#endif // GRAPH_H
