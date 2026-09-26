"""Keep the sidebars from repeating a parent's name.

Two jobs, one rule: an entry says its own name and not its parent's.
``rewrite`` nests the left sidebar by dotted name; ``shorten`` strips the
parent prefix from the right-hand "On this page" panel.

Nesting the left sidebar
------------------------

``autosummary`` writes one flat ``toctree`` entry per documented object, so
the sidebar lists ``rotifer.db.ncbi.entrez`` as a sibling of ``rotifer.db``
instead of a child of it, and repeats the parent name in every label. This
extension rewrites the global table of contents after Sphinx renders it:
sibling entries whose titles are dotted Python names are regrouped into the
tree those names describe, and each label is shortened to its own last
component.

The rewrite is HTML post-processing only. No source file, ``toctree`` or
``autosummary`` block changes, the ``localtoc`` on the right is untouched,
and nothing here runs for a non-HTML builder.

A package that has no page of its own still gets a node: ``rotifer.core.io``
is reachable only through ``rotifer.core.io.ftp``, so it becomes a label
that toggles its branch rather than a link. Shibuya's sidebar script binds
its collapse handler to any anchor whose ``href`` is ``#``, which is also
how the theme renders the current page, so those labels fold and unfold
like every other node.

Collapsing itself needs no help here. Shibuya appends a toggle button to
every sidebar item that has children, so deepening the tree is what gives
each module and each page its own toggle.

Shortening the "On this page" panel
-----------------------------------
autodoc labels a method entry with the class it belongs to, so a class page
listed ``NucleotideFeaturesCursor`` and then, indented under it,
``NucleotideFeaturesCursor.fetchall()``. Those labels are one unbreakable
token wider than the panel, so they were clipped mid-word against its right
edge. ``shorten`` drops the prefix a child already inherits from its
parent, leaving ``fetchall()``, which fits. The stylesheet still allows a
long label to wrap, so nothing is ever clipped again.
"""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from typing import Any, Dict, List, Optional

# A title that is safe to split on dots: a Python dotted name and nothing
# else. Section titles ("API reference", "Data access") never match, which
# is what keeps the rewrite inside the generated part of the tree.
DOTTED_NAME = re.compile(r"^[A-Za-z_]\w*(?:\.\w+)*$")

# Marks the two kinds of node the stylesheet needs to tell apart: a label
# that only toggles its branch, and an item with no children, which needs no
# room for a toggle button.
BRANCH_CLASS = "rot-toc-branch"
LEAF_CLASS = "rot-toc-leaf"


class _Node:
    """One sidebar item: a label, the page it links to, and its children."""

    __slots__ = ("label", "anchor", "children", "current")

    def __init__(self, label: str) -> None:
        self.label = label
        self.anchor: Optional[ET.Element] = None
        self.children: Dict[str, "_Node"] = {}
        self.current = False


def _text_of(element: ET.Element) -> str:
    return "".join(element.itertext()).strip()


def _classes_of(element: ET.Element) -> List[str]:
    return (element.get("class") or "").split()


def _is_dotted_group(ul: ET.Element) -> bool:
    """True when every item of ``ul`` is a dotted name and at least one nests."""
    items = ul.findall("li")
    if not items:
        return False
    nested = False
    for li in items:
        anchor = li.find("a")
        if anchor is None:
            return False
        title = _text_of(anchor)
        if not DOTTED_NAME.match(title):
            return False
        nested = nested or "." in title
    return nested


def _collect(ul: ET.Element, found: List[ET.Element]) -> None:
    """Flatten a sidebar branch into its anchors, at any depth."""
    for li in ul.findall("li"):
        anchor = li.find("a")
        if anchor is not None:
            found.append(anchor)
        child = li.find("ul")
        if child is not None:
            _collect(child, found)


def _build(anchors: List[ET.Element]) -> _Node:
    root = _Node("")
    for anchor in anchors:
        current = "current" in _classes_of(anchor)
        node = root
        for part in _text_of(anchor).split("."):
            node = node.children.setdefault(part, _Node(part))
            # Mark the whole path so the theme opens the tree down to the
            # page being read.
            node.current = node.current or current
        # The same name can arrive twice, once as a module listed in
        # docs/api/index.rst and once as a member of its parent's page. The
        # first anchor wins; both point at the same file.
        if node.anchor is None:
            node.anchor = anchor
    return root


def _order(item):
    """Packages first, then pages, each group alphabetical and case-blind.

    The curated order of the autosummary blocks in docs/api/index.rst cannot
    survive the regroup, because it is an order over a flat list and this is
    a tree. Sorting instead puts every name in the same place on every page.
    """
    label, node = item
    return (0 if node.children else 1, label.lower(), label)


def _anchor_for(node: _Node) -> ET.Element:
    anchor = ET.Element("a")
    if node.anchor is not None:
        for name, value in node.anchor.attrib.items():
            anchor.set(name, value)
    else:
        # A package with no page of its own. "#" is what makes the theme
        # treat the label as a toggle instead of a link.
        anchor.set("class", "reference internal " + BRANCH_CLASS)
        anchor.set("href", "#")
    anchor.text = node.label
    return anchor


def _emit(node: _Node, level: int, into: ET.Element) -> None:
    for _, child in sorted(node.children.items(), key=_order):
        li = ET.SubElement(into, "li")
        classes = ["toctree-l%d" % level]
        if child.current:
            classes.append("current")
        if not child.children:
            classes.append(LEAF_CLASS)
        li.set("class", " ".join(classes))
        li.append(_anchor_for(child))
        if child.children:
            _emit(child, level + 1, ET.SubElement(li, "ul"))


def _regroup(ul: ET.Element, level: int) -> None:
    anchors: List[ET.Element] = []
    _collect(ul, anchors)
    tree = _build(anchors)
    for child in list(ul):
        ul.remove(child)
    ul.text = None
    _emit(tree, level, ul)


def _walk(ul: ET.Element, level: int) -> None:
    if _is_dotted_group(ul):
        _regroup(ul, level)
        return
    for li in ul.findall("li"):
        child = li.find("ul")
        if child is not None:
            _walk(child, level + 1)


def rewrite(toc: str) -> str:
    """Return ``toc`` with its dotted-name branches nested and relabelled."""
    if not toc or not toc.strip():
        return toc
    try:
        root = ET.fromstring("<div>" + toc.strip() + "</div>")
    except ET.ParseError:
        # Never break a build over the sidebar: fall back to the flat tree.
        return toc
    lists = root.findall("ul")
    if not lists:
        return toc
    for ul in lists:
        _walk(ul, 1)
    return ET.tostring(root, encoding="unicode")[len("<div>"):-len("</div>")]


def _sole_text_node(element: ET.Element) -> Optional[ET.Element]:
    """The one descendant that holds all of ``element``'s text, if there is one.

    An "On this page" label is normally ``a > code > span`` with the text in
    the span. Anything else, such as a label split across inline markup, is
    left alone rather than guessed at.
    """
    holders = [node for node in element.iter() if (node.text or "").strip()]
    if len(holders) != 1:
        return None
    return holders[0]


def _shorten_children(ul: ET.Element, parent: str) -> None:
    for li in ul.findall("li"):
        anchor = li.find("a")
        label = _text_of(anchor) if anchor is not None else ""
        if anchor is not None and parent and " " not in parent:
            if label.startswith(parent + "."):
                holder = _sole_text_node(anchor)
                if holder is not None:
                    holder.text = holder.text.strip()[len(parent) + 1:]
                    label = _text_of(anchor)
        child = li.find("ul")
        if child is not None:
            # Match against the label as autodoc wrote it, which is what the
            # grandchildren repeat, not against the shortened one.
            _shorten_children(child, _text_of(anchor) if anchor is None else label)


def shorten(toc: str) -> str:
    """Return ``toc`` with each entry's inherited parent prefix removed."""
    if not toc or not toc.strip():
        return toc
    try:
        root = ET.fromstring("<div>" + toc.strip() + "</div>")
    except ET.ParseError:
        return toc
    lists = root.findall("ul")
    if not lists:
        return toc
    for ul in lists:
        _shorten_children(ul, "")
    return ET.tostring(root, encoding="unicode")[len("<div>"):-len("</div>")]


def _patch_page_context(
    app: Any,
    pagename: str,
    templatename: str,
    context: Dict[str, Any],
    doctree: Any,
) -> None:
    """Rewrite the left sidebar's callable and the right panel's markup."""
    render = context.get("toctree")
    if callable(render):

        def render_nested(**kwargs: Any) -> str:
            return rewrite(render(**kwargs))

        context["toctree"] = render_nested

    toc = context.get("toc")
    if isinstance(toc, str):
        context["toc"] = shorten(toc)


def setup(app: Any) -> Dict[str, Any]:
    app.connect("html-page-context", _patch_page_context)
    return {
        "version": "1.0",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
