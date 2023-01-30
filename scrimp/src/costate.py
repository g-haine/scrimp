from ..src.state import State


class CoState(State):
    def __init__(self, name: str, description: str, state: State, substituted=False):
        assert isinstance(state, State), (f'State {state} must be added before its co-state')
        super().__init__(name, description, state.get_kind(), state.get_region(), state.get_mesh_id())
        self._state = state
        self._substituted = substituted

    def __str__(self):
        where = ''
        if self.get_state().get_region() is not None:
            where = ', in region numbered ' + str(self.get_state().get_region())
        print('A co-state variable', self.get_name(), ', describing \'', self.get_description(),
              '\', associated to state \'', self.get_state(),
              '\' has been initialized as a', self.get_kind(), 'on mesh', self.get_mesh_id(),
              where)
        print('The constitutive relations between the state', self.get_state(), 'and the co-state', self.get_name(),
              'will' + (not self.get_substituted()) * ' not' + ' be substituted for the resolution: variable',
              self.get_name(),
              'will' + self.get_substituted() * ' not' + ' be considered as an unknown')

    def get_state(self):
        return self._state

    def get_substituted(self):
        return self._substituted

    #
    # def add_costate(self, name, description, state, substituted=False):
    #     """
    #     Add a variable to the dict `costates` of the dpHs
    #
    #     A `state` is mandatory. The kind of field of the co-state is that of the state
    #
    #     A `port` named as the associated state is created
    #
    #     :param name: the name of the costate
    #     :type name: str
    #     :param description: a physically motivated description (e.g. `velocity`)
    #     :type description: str
    #     :param state: the name of the state
    #     :type state: str
    #     :param substituted: if 'True' (default: `False`) the constitutive relations are substituted into the dynamic
    #     :type substituted: bool
    #
    #     :return:
    #         * append the new costate to the dict `costates` of the dpHs
    #         * append the non-algebraic port (state, costate) to the dict `ports` of the dpHs
    #     """
    #
    #
    #
    #     self.costates[name] = {'description': description,
    #                            'kind': self.states[state]['kind'],
    #                            'region': self.states[state]['region'],
    #                            'mesh_id': self.states[state]['mesh_id'],
    #                            'state': state,
    #                            'port': None
    #                            }
    #
    #     self.states[state]['costate'] = name
    #
    #     self.add_port(state, state, name, self.costates[name]['kind'],
    #                   self.states[state]['mesh_id'], algebraic=False,
    #                   substituted=substituted,
    #                   region=self.states[state]['region'])
    #
    #     self.states[state]['port'] = state
    #     self.costates[name]['port'] = state
    #
    #     where = ''
    #     if self.states[state]['region'] is not None:
    #         where = ', in region numbered ' + str(self.states[state]['region'])
    #     print('A co-state variable', name, ', describing \'', description, '\', associated to state \'', state,
    #           '\' has been initialized as a', self.costates[name]['kind'], 'on mesh', self.costates[name]['mesh_id'],
    #           where)
    #     print('The constitutive relations between the state', state, 'and the co-state', name,
    #           'will' + (not substituted) * ' not' + ' be substituted for the resolution: variable', name,
    #           'will' + substituted * ' not' + ' be considered as an unknown')
