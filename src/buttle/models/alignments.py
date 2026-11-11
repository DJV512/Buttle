from dataclasses import dataclass
from . import molecules as btmm
from collections import Counter

@dataclass
class AlignedPair:
    sequence1: btmm.Protein | btmm.DNA
    sequence2: btmm.Protein | btmm.DNA
    align1: list[str]
    align2: list[str]
    distance: float


@dataclass
class Cluster:
    id: str
    name: str
    left: btmm.Protein | btmm.DNA | Cluster
    right: btmm.Protein | btmm.DNA | Cluster
    members: list[btmm.Protein] | list[btmm.DNA]

    def __str__(self):
        members = ", ".join(member.name for member in self.members)
        return f"Cluster(name={self.name}, left={self.left}, right={self.right}, members={members})"


    def swap_nodes(self) -> None:
        left_node = self.left
        right_node = self.right
        self.left = right_node
        self.right = left_node

    def print_tree(self) -> None:
        """
        Print the guide tree using Unicode box-drawing characters.
        """
        print("root")
        self._print_children("")


    def _print_children(self, prefix: str) -> None:
        """
        Recursively print the left and right children of a cluster.
        """
        for child, is_last in ((self.left, False), (self.right, True)):
            connector = "└── " if is_last else "├── "
            label = "*" if isinstance(child, Cluster) else child.name

            print(f"{prefix}{connector}{label}")

            if isinstance(child, Cluster):
                extension = "    " if is_last else "│   "
                child._print_children(prefix + extension)


    
@dataclass
class Profile:
    name: list[str]
    alignment: list[list[str]]
    molecule_type: str 

    def __str__(self):
        text = ""
        for name, align in zip(self.name, self.alignment):
            text += f"{name}, {align}\n"

        return text

    def add_consensus(self):
        consensus = []
        
        aligned_length = len(self.alignment[0])
    
        for i in range(aligned_length):
            residues = []
            for alignment in self.alignment:
                residues.append(alignment[i])
            count = Counter(residues)
            consensus.append(count.most_common(1)[0][0])
    
        self.name.append("Consensus")
        self.alignment.append(consensus)



@dataclass
class ParsedProfile:
    name: list[str]
    alignment: list[list[tuple[str, str]]]