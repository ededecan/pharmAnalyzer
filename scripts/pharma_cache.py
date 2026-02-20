"""
Caching and indexing utilities for pharmacogenomics analyzer.
Improves performance through memoization and pre-computed indices.
"""

from typing import Dict, List, Set, Tuple, Any, Callable
from functools import lru_cache, wraps
import time


class CacheManager:
    """Manages caching for expensive operations"""
    
    def __init__(self):
        self.cache = {}
        self.hits = 0
        self.misses = 0
    
    def get(self, key: str, compute_fn: Callable, *args, **kwargs) -> Any:
        """Get cached value or compute and cache it"""
        if key in self.cache:
            self.hits += 1
            return self.cache[key]
        
        self.misses += 1
        value = compute_fn(*args, **kwargs)
        self.cache[key] = value
        return value
    
    def clear(self):
        """Clear cache"""
        self.cache.clear()
        self.hits = 0
        self.misses = 0
    
    def stats(self) -> Dict[str, int]:
        """Get cache statistics"""
        total = self.hits + self.misses
        hit_rate = (self.hits / total * 100) if total > 0 else 0
        return {
            'hits': self.hits,
            'misses': self.misses,
            'total': total,
            'hit_rate': hit_rate
        }


class RsIDIndex:
    """Pre-computed index for rapid rsID lookups"""
    
    def __init__(self, allele_definitions: Dict[str, Dict[str, Any]]):
        """Build index from allele definitions
        
        Index structure:
        {
            'rsid': {
                'gene': {
                    'allele': [nucleotide_values]
                }
            }
        }
        """
        self.index = {}
        self._build_index(allele_definitions)
    
    def _build_index(self, allele_definitions: Dict[str, Dict[str, Any]]):
        """Build inverted index from rsID to alleles"""
        for gene, gene_data in allele_definitions.items():
            alleles = gene_data.get('alleles', {})
            
            for allele_name, allele_values in alleles.items():
                for rsid in allele_values.keys():
                    if rsid not in self.index:
                        self.index[rsid] = {}
                    if gene not in self.index[rsid]:
                        self.index[rsid][gene] = {}
                    
                    # Store allele name that has this rsID
                    self.index[rsid][gene][allele_name] = allele_values[rsid]
    
    def lookup_gene_alleles(self, rsid: str, gene: str) -> List[str]:
        """Get all alleles for a gene that have this rsID"""
        if rsid not in self.index:
            return []
        if gene not in self.index[rsid]:
            return []
        return list(self.index[rsid][gene].keys())
    
    def has_rsid(self, rsid: str) -> bool:
        """Check if rsID exists in index"""
        return rsid in self.index


class DiplotypeLookupCache:
    """Pre-computed cache for diplotype lookups"""
    
    def __init__(self, diplotype_phenotypes: Dict[str, Dict[str, Dict[str, str]]]):
        """Build cache from diplotype phenotypes
        
        Cache structure:
        {
            'gene': {
                'diplotype_exact': {...phenotype_info...},
                'diplotype_reversed': {...phenotype_info...}
            }
        }
        """
        self.cache = {}
        self._build_cache(diplotype_phenotypes)
    
    def _build_cache(self, diplotype_phenotypes: Dict[str, Dict[str, Dict[str, str]]]):
        """Build bidirectional diplotype cache"""
        for gene, diplotypes in diplotype_phenotypes.items():
            self.cache[gene] = {}
            
            for diplotype, phenotype_info in diplotypes.items():
                # Store original
                self.cache[gene][diplotype] = phenotype_info
                
                # Store reversed variant
                parts = diplotype.split('/')
                if len(parts) == 2:
                    reversed_dipl = f"{parts[1]}/{parts[0]}"
                    if reversed_dipl != diplotype:
                        self.cache[gene][reversed_dipl] = phenotype_info
    
    def lookup(self, gene: str, diplotype: str) -> Dict[str, str]:
        """Lookup diplotype with single O(1) operation"""
        if gene not in self.cache:
            return {}
        
        return self.cache[gene].get(diplotype, {})


class AlleleFrequencyIndex:
    """Quick lookup for allele frequencies"""
    
    def __init__(self, allele_frequencies: Dict[str, Dict[str, str]]):
        """Index allele frequencies by gene and allele"""
        self.index = allele_frequencies
    
    def get_frequency(self, gene: str, allele: str, default: str = "N/A") -> str:
        """Get frequency in O(1) time"""
        return self.index.get(gene, {}).get(allele, default)


def time_operation(func: Callable) -> Callable:
    """Decorator to time function execution"""
    @wraps(func)
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        elapsed = time.time() - start
        return result, elapsed
    return wrapper


class VariantGrouper:
    """Efficient grouping of variants by gene"""
    
    @staticmethod
    def group_by_gene(variants: List[Dict[str, str]]) -> Dict[str, List[Dict[str, str]]]:
        """Group variants by gene in single pass"""
        groups = {}
        for variant in variants:
            gene = variant.get('symbol', variant.get('gene', 'Unknown'))
            if gene not in groups:
                groups[gene] = []
            groups[gene].append(variant)
        return groups
    
    @staticmethod
    def group_by_rsid(variants: List[Dict[str, str]]) -> Dict[str, List[Dict[str, str]]]:
        """Group variants by rsID for quick lookup"""
        groups = {}
        for variant in variants:
            rsid = variant.get('existing_variation', '-')
            if rsid and rsid != '-':
                if rsid not in groups:
                    groups[rsid] = []
                groups[rsid].append(variant)
        return groups
