import time
import numpy as np
import matplotlib.pyplot as plt
import getfem as gf

from src.domain import Domain


class Term:
    """This class defines a term for the Hamiltoninan."""

    def __init__(
        self, description: str, expression: str, regions: str, mesh_id: int = 0
    ):
        """This constructor defines the object Term for the Hamiltonian functional.

        Args:
            description (str): the name or description of the term (e.g. 'Kinetic energy')
            expression (str): the formula, using the `Model` variables, defining the term. Parameters are allowed (e.g. '0.5*q.T.q')
            regions (str): the region ids of the mesh mesh_id where the expression has to be evaluated
            mesh_id (sint): the mesh id of the mesh where the regions belong
        """
        self.__description = description
        self.__expression = expression
        self.__regions = regions
        self.__mesh_id = mesh_id
        self.__values = []

    def get_description(self) -> str:
        """This function gets the description of the term.

        Returns:
            str: description of the term
        """
        return self.__description

    def get_expression(self) -> str:
        """This function gets the matematical expression of the term.

        Returns:
            str: matematical expression of the term.
        """
        return self.__expression

    def get_regions(self) -> str:
        """This function gets the regions of the term.

        Returns:
            str: regions of the termthe region ids of the mesh mesh_id where the expression has to be evaluated
        """
        return self.__regions

    def get_mesh_id(self) -> int:
        """This function gets the mesh id of the mesh where the regions belong.

        Returns:
            int: the mesh id of the mesh where the regions belong
        """
        return self.__mesh_id

    def get_values(self) -> list:
        """This function gets the valeus of the term.

        Returns:
            list: list of values of the term
        """
        return self.__values.copy()

    def set_value(self, value):
        """This function sets a valeus for the term.

        Args:
            value : a value for the term
        """

        # TODO assert to check the type of the vlaue
        self.__values.append(value)


class Hamiltonian:
    def __init__(self, name: str) -> None:
        self.__terms = []
        self.__name = name
        self.__is_computed = False

    def add_term(self, term: Term):
        """This function add a term to the term list of the Hamiltonian

        Args:
            term (Term): term for the Hamiltonian
        """
        assert isinstance(term, Term)
        self.__terms.append(term)

    def compute(self, domain: Domain, solution: dict):
        """
        Compute each `term` constituting the Hamiltonian

        :return: fill the Hamiltonian[term]['values'] with a list of values at times t
        """

        assert (
            self.solve_done
        ), "System has not been solved yet, Hamiltonian can not be computed"

        print("Start computing the Hamiltonian")
        start = time.perf_counter()
        for t in range(len(solution["t"])):
            self.gf_model.to_variables(solution["z"][t])
            for term in self.__terms:
                term_value_at_t = 0.0
                for region in term.get_regions():
                    term_value_at_t += gf.asm(
                        "generic",
                        domain.int_method[term.get_mesh_id()],
                        0,
                        term.get_expression(),
                        region,
                        self.gf_model,
                    )
                self.add_term(term_value_at_t)

        print("Hamiltonian has been computed in", time.perf_counter() - start, "s")
        self.__is_computed = True

    def plot_Hamiltonian(self, solution: dict, with_powers=True):
        """
        Plot each term constituting the Hamiltonian and the Hamiltonian

        May include the power terms

        :param with_powers: if `True` (default), the plot wil also contains the power of each algebraic ports
        :type with_powers: bool

        :return: a matplotlib figure
        """

        if not self.__is_computed:
            self.compute()

        t = np.array(solution["t"])
        fig = plt.figure(figsize=[8, 5])
        ax = fig.add_subplot(111)
        HamTot = np.zeros(t.size)

        for term in self.__terms:
            values = np.array(term.get_values())
            HamTot += values
            ax.plot(t, values, label=term.get_description())

        if len(self.__terms) > 1:
            ax.plot(t, HamTot, label=self.__name)

        if with_powers:
            self.plot_powers(ax, HamTot=HamTot)

        ax.legend()
        ax.grid(axis="both")
        ax.set_xlabel("time t")
        ax.set_ylabel("Hamiltonian terms")
        ax.set_title("Evolution of Hamiltonian terms")
        plt.show()

    def get_name(self) -> str:
        """This function returns the name of the Hamiltonion.

        Args:
            name (str): name of the Hamiltonian
        """

        return self.__name

    def set_name(self, name: str):
        """This function set the name for the Hamiltonion. For plotting purposes.

        Args:
            name (str): name of the Hamiltonian
        """

        self.__name = name

    def __len__(self) -> int:
        """This function returns the length of the Hamiltonian, that is the number of therms.

        Returns:
            int: _description_
        """
        return len(self.__terms)

    def get_terms(self) -> list:
        """This function returns a copy of the list of terms of the Hamiltonian.

        Returns:
            list: list of terms of the Hamiltonian.
        """
        return self.__terms.copy()

    def __getitem__(self, index) -> Term:
        """This function return the term correspondig ti the index.

        Args:
            index (_type_): _description_

        Returns:
            Term: _description_
        """
        try:
            return self.__terms[index]
        except IndexError:
            print(f"index {index} out of {len(self.__terms)}")

    def set_is_computed(self):
        """This function sets the Hamiltonian as computed."""
        self.__is_computed = True

    def get_is_computed(self) -> bool:
        """This function returns True if the Hamiltonian is computed, False otherwise.

        Returns:
            bool: _description_
        """
        return self.__is_computed
