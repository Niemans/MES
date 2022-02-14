# MES
Obliczanie temperatury w każdym węźle prostokątnej siatki BxH w podanym przedziale czasowym.<br />
Plik MES.cpp jest głównym plikiem z funkcją main.<br />

zmienne ogólne:<br />
gaussNumber - liczba punktów całkowania gaussa potrzebna do interpolacji całkowania metodą gaussa (działa tylko dla 2 i 3 punktów, gdyż dla innych liczb nie ma wpisanych wartości)<br />
dimentions - liczba wymiarów obliczania siatki, lecz działa tylko przy zmiennej równej 2<br />

stałe danego obiektu:<br />
k - współczynnik przewodzenia ciepła<br />
a - współczynnik konwekcji (alfa)<br />
c - pojemność cieplna<br />
ro - gęstość<br />


zmienne czasu:<br />
tau - krok czasowy<br />
tend - czas końcowy<br />

zmienne temperatury:<br />
T0 - temperatura początkowa obiektu<br />
Tot - temperaturaotoczenia<br />

zmienne określający siatkę:<br />
b - szerokość siatki<br />
h - wysokość siatki<br />
nB - liczba węzłów po szerokości<br />
nH - liczba węzłów po wysokości<br />

Każda funkcja w plikach z prefixem "print" można odkomentować i sprawdzić, jakie wyniki wychodzą w danym momencie działania programu.<br />
Reszta może być mniej bezpieczna do zmiany.
