{{ fullname | escape | underline }}

.. automodule:: {{ fullname }}

.. currentmodule:: {{ fullname }}

{% set visible_functions = functions | reject("in", skip_members) | list %}
{% set visible_classes = classes | reject("in", skip_members) | list %}
{% set visible_exceptions = exceptions | reject("in", skip_members) | list %}

{% block functions %}
{% if visible_functions %}
.. rubric:: Functions

.. autosummary::
   :toctree:
{% for item in visible_functions %}
   {{ item }}
{%- endfor %}
{% endif %}
{% endblock %}

{% block classes %}
{% if visible_classes %}
.. rubric:: Classes

.. autosummary::
   :toctree:
{% for item in visible_classes %}
   {{ item }}
{%- endfor %}
{% endif %}
{% endblock %}

{% block exceptions %}
{% if visible_exceptions %}
.. rubric:: Exceptions

.. autosummary::
   :toctree:
{% for item in visible_exceptions %}
   {{ item }}
{%- endfor %}
{% endif %}
{% endblock %}
